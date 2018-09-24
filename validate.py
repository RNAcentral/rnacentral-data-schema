#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import os
import json
import logging
from collections import Counter

import click
import requests
import jsonschema as js
# import jsonref as jr


HERE = os.path.abspath(os.path.dirname(__file__))
SECTIONS = os.path.join(HERE, 'sections')
SCHEMA_NAME = 'rnacentral-schema.json'
LOGGER = logging.getLogger(__name__)

TAX_URL = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{taxon_id}'

PUB_URLS = {
    'PMID': 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{value}+AND+SRC:MED&format=json',
    'DOI': 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:{value}&format=json',
}


class ValidationWarning(js.ValidationError):
    pass


class UnrunnableValidator(js.ValidationError):
    pass


class ExtendedValidator(js.validators.Draft4Validator):
    def __init__(self, schema, extra, *args, **kwargs):
        self.extra_validators = extra
        super(ExtendedValidator, self).__init__(schema, *args, **kwargs)

    def update_error(self, instance, validator, item, error):
        # Copied from the jsonschema validators
        try:
            name = validator.__name__
        except AttributeError:
            name = validator.__class__.__name__

        error._set(
            validator=name,
            validator_value=item,
            instance=instance,
            schema=None
        )
        return error

    def validate_metadata(self, instance):
        for validator in self.extra_validators:
            if not hasattr(validator, 'validate_metadata'):
                continue
            try:
                value = instance['metaData']
                for error in validator.validate_metadata(value):
                    yield (validator, value, error)
            except Exception as err:
                yield (validator, value, UnrunnableValidator(err))

    def validate_ncrnas(self, instance):
        for ncrna in instance['data']:
            for validator in self.extra_validators:
                if not hasattr(validator, 'validate_ncrna'):
                    continue
                try:
                    for error in validator.validate_ncrna(ncrna):
                        yield (validator, ncrna, error)
                except Exception as err:
                    yield (validator, ncrna, UnrunnableValidator(err))

    def iter_errors(self, instance, _schema=None):
        parent = super(ExtendedValidator, self)
        for error in parent.iter_errors(instance, _schema=_schema):
            yield error

        if isinstance(instance, dict):
            if 'metaData' in instance:
                for validator, value, error in self.validate_metadata(instance):
                    yield self.update_error(instance, validator, value, error)

            if 'data' in instance:
                for validator, value, error in self.validate_ncrnas(instance):
                    yield self.update_error(instance, validator, value, error)


class KnownGlobalIdValidator(object):
    def __init__(self):
        with open('sections/data-provider.json', 'r') as raw:
            known = json.load(raw)
            self.known = set(known["properties"]['dataProvider']['enum'])

    def validate_ncrna(self, ncrna):
        gene_id = ncrna.get('gene', {}).get('geneId', None)
        if gene_id:
            name, _ = gene_id.split(':', 1)
            if name.upper() not in self.known:
                yield js.ValidationError("Unknown database: %s" % name)

        for global_id in ncrna.get('crossReferenceIds', []):
            name, _ = global_id.split(':', 1)
            if name.upper() not in self.known:
                yield js.ValidationError("Xref to unknown db: %s" % name)


class ActiveTaxonIdValidator(object):
    def __init__(self):
        self.seen = set()
        self.failed = set()

    def validate_ncrna(self, ncrna):
        taxon_id = int(ncrna['taxonId'].split(':', 1)[1])
        if taxon_id in self.seen:
            return

        if taxon_id in self.failed:
            yield js.ValidationError("Invalid Taxon id: %s" % taxon_id)
        else:
            try:
                response = requests.get(TAX_URL.format(taxon_id=taxon_id))
                response.raise_for_status()
                self.seen.add(taxon_id)
            except requests.HTTPError:
                self.failed.add(taxon_id)
                yield js.ValidationError("Invalid Taxon id: %s" % taxon_id)


class PublicationValidator(object):
    def __init__(self):
        self.seen = set()
        self.failed = set()
        self.requires_ncrna_publications = True

    def validate_pmid(self, pub_id):
        db, db_id = pub_id.split(':', 1)
        if db_id in self.seen:
            return

        if db not in PUB_URLS:
            return ValidationWarning("Could not validate %s" % pub_id)

        if db_id in self.failed:
            return js.ValidationError("Invalid reference id: %s" % pub_id)
        try:
            url = PUB_URLS[db].format(value=db_id)
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            if data['hitCount'] == 0:
                raise requests.HTTPError("Not found")
            self.seen.add(db_id)
        except requests.HTTPError:
            self.failed.add(db_id)
            return js.ValidationError("Invalid reference id: %s" % pub_id)

    def validate_metadata(self, metadata):
        publications = metadata.get('publications', [])
        if not publications:
            yield ValidationWarning("Databases should have a reference")
        else:
            self.requires_ncrna_publications = False

        for pmid in publications:
            error = self.validate_pmid(pmid)
            if error:
                yield error

    def validate_ncrna(self, ncrna):
        publications = ncrna.get('publications', [])
        if not publications and self.requires_ncrna_publications:
            yield js.ValidationError("Must have at least one reference")

        for pub_id in publications:
            error = self.validate_pmid(pub_id)
            if error:
                yield error


class SecondaryStructureValidator(object):
    def validate_ncrna(self, ncrna):
        if 'secondaryStructure' not in ncrna:
            return

        if len(ncrna['secondaryStructure']) != len(ncrna['sequence']):
            yield js.ValidationError("Secondary structure wrong size")


# Not clear if this is actually needed.
# def trna_annotations(ncrna):
#     isoType = ncrna.get('additionalAnnotations', {}).get('isoType', None)
#     anticodon = ncrna.get('sequenceFeatures', {}).get('anticodon', None)
#     if isoType or anticodon:
#         if ncrna['soTermId'] != 'SO:0000253':
#             yield js.ValidationError("tRNA has the wrong SO term")


class NameValidator(object):
    def validate_ncrna(self, ncrna):
        name = None
        if 'description' in ncrna and ncrna['description']:
            LOGGER.debug("Using transcript description for name of %s",
                         ncrna['primaryId'])
            name = ncrna['description']
        if 'name' in ncrna and ncrna['name']:
            LOGGER.debug("Using transcript name for name of %s",
                         ncrna['primaryId'])
            name = ncrna['name']
        if 'gene' in ncrna:
            gene = ncrna['gene']
            if 'name' in gene:
                LOGGER.debug("Using gene name for name of %s",
                             ncrna['primaryId'])
                name = gene['name']
            if 'symbol' in gene:
                LOGGER.debug("Using gene symbol for name of %s",
                             ncrna['primaryId'])
                name = gene['symbol']

        if name:
            LOGGER.debug("Using name %s for %s", name, ncrna['primaryId'])
        else:
            yield js.ValidationError("No name for %s" % ncrna['primaryId'])


class CoordinateDirectionValidator(object):
    def validate_ncrna(self, ncrna):
        for location in ncrna.get('genomeLocations', []):
            for exon in location['exons']:
                if exon['strand'] == '+' or exon['strand'] == '.' or \
                        exon['strand'] == '-':
                    if not exon['startPosition'] < exon['endPosition']:
                        yield js.ValidationError("Start must be < end: %s")
                else:
                    raise ValueError("Shouldn't be here")


class AcceptableUncertaintyValidator(object):
    def validate_ncrna(self, ncrna):
        standard = set('ACGT')
        sequence = ncrna['sequence']
        total = float(len(ncrna['sequence']))
        uncertainty = sum(1 for s in sequence if s not in standard)
        if float(uncertainty) / total > 0.1:
            yield ValidationWarning("Sequence for %s is too uncertain (%f/%i)" %
                                    (ncrna, uncertainty, total))


def validate(data, schema_path, sections_path):

    with open(schema_path, 'r') as raw:
        schema = json.load(raw)

    validators = [
        AcceptableUncertaintyValidator(),
        CoordinateDirectionValidator(),
        NameValidator(),
        SecondaryStructureValidator(),
        ActiveTaxonIdValidator(),
        KnownGlobalIdValidator(),
        PublicationValidator(),
    ]

    base = 'file://%s/' % sections_path
    validator = ExtendedValidator(
        schema,
        validators,
        format_checker=js.FormatChecker(),
        resolver=js.RefResolver(base, None),
    )

    found = False
    counts = Counter()
    for error in validator.iter_errors(data):
        counts[error.validator] += 1
        if isinstance(error, ValidationWarning):
            LOGGER.warning(error.message)
        else:
            found = True
            LOGGER.error(error.message)

    if found:
        summary = ', '.join('%s: %s' % (k, v) for k, v in counts.items())
        raise click.ClickException("Validation failed: %s" % summary)


@click.command()
@click.argument('filename')
@click.option('--schema', default=SCHEMA_NAME,
              help='Filename of the schema to use')
@click.option('--sections', default=SECTIONS,
              help='Directory where schema parts are kept')
def main(filename, schema=None, sections=None):
    with open(filename, 'r') as raw:
        data = json.load(raw)

    data['metaData']["dataProvider"] = data['metaData']["dataProvider"].upper()
    validate(data, schema, os.path.abspath(sections))


if __name__ == '__main__':

    logging.basicConfig(
        format='%(levelname)s: %(message)s',
        level=logging.WARNING,
    )

    main()
