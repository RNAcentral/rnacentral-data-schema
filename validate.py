#!/usr/bin/env python

import os
import json
import logging

import click
import jsonschema as js
# import jsonref as jr


HERE = os.path.abspath(os.path.dirname(__file__))
SECTIONS = os.path.join(HERE, 'sections')
SCHEMA_NAME = 'rnacentral-schema.json'
LOGGER = logging.getLogger(__name__)


def validate_secondary_structure(data):
    for ncrna in data['data']:
        if 'secondaryStructure' not in ncrna:
            continue
        assert(len(ncrna['secondaryStructure']) == len(ncrna['sequence']))


def validate_is_known_global_ids(data):
    with open('sections/data-provider.json', 'r') as raw:
        known = json.load(raw)
        known = set(known["properties"]['dataProvider']['enum'])

    for ncrna in data["data"]:
        gene_id = ncrna.get('gene', {}).get('geneId', None)
        if gene_id:
            name, _ = gene_id.split(':', 1)
            assert name in known, "Unknown database: %s" % name

        for global_id in ncrna.get('crossReferenceIds', []):
            name, _ = global_id.split(':', 1)
            assert name in known, "Xref to unknown db: %s" % name


# Not clear if this is actually needed.
def validate_trna_annotations(data):
    for ncrna in data['data']:
        isoType = ncrna.get('additionalAnnotations', {}).get('isoType', None)
        anticodon = ncrna.get('sequenceFeatures', {}).get('anticodon', None)
        if isoType or anticodon:
            assert ncrna['soTermId'] == 'SO:0000253'


# Unsure if we should require this, maybe make it an option? I will leave the
# code here for now.
def validate_id_format(data):
    expected = data['metaData']['dataProvider']
    for ncrna in data['data']:
        primary_id = ncrna['primaryId']
        db, _ = primary_id.split(':', 1)
        if db != expected:
            msg = "Expected %s to start with %s" % (primary_id, expected)
            raise js.ValidationError(msg)

        gene_id = ncrna.get('gene', {}).get('geneId', None)
        if gene_id:
            gene_db = primary_id.split(':', 1)
            assert gene_db == expected


def validate_can_produce_name(data):
    for ncrna in data['data']:
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
            raise ValueError("No name for %s", ncrna['primaryId'])


def validate_coordinate_direction(data):
    for ncrna in data['data']:
        for location in ncrna.get('genomeLocations', []):
            for exon in location['exons']:
                if exon['strand'] == '+' or exon['strand'] == '.' or \
                        exon['strand'] == '-':
                    assert exon['startPosition'] < exon['endPosition']
                else:
                    raise ValueError("Shouldn't be here")


def validate(data, schema_path, sections_path):

    with open(schema_path, 'r') as raw:
        schema = json.load(raw)

    base = 'file://%s/' % sections_path
    js.validate(
        data,
        schema,
        format_checker=js.FormatChecker(),
        resolver=js.RefResolver(base, None),
    )
    validate_secondary_structure(data)
    validate_is_known_global_ids(data)
    validate_trna_annotations(data)
    validate_can_produce_name(data)
    validate_coordinate_direction(data)


@click.command()
@click.argument('filename')
@click.option('--schema', default=SCHEMA_NAME,
              help='Filename of the schema to use')
@click.option('--sections', default=SECTIONS,
              help='Directory where schema parts are kept')
def main(filename, schema=None, sections=None):
    with open(filename, 'r') as raw:
        data = json.load(raw)

    validate(data, schema, os.path.abspath(sections))


if __name__ == '__main__':
    main()
