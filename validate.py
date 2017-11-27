#!/usr/bin/env python

import os
import json

import click
import jsonschema as js


HERE = os.path.abspath(os.path.dirname(__file__))
SECTIONS = os.path.join(HERE, 'sections')


def validate_secondary_structure(data):
    for entry in data['data']:
        if 'secondaryStructure' not in entry:
            continue
        assert(len(entry['secondaryStructure']) == len(entry['sequence']))


def validate_is_known_global_ids(data):
    with open('sections/data-provider.json', 'rb') as raw:
        known = json.load(raw)
        known = set(known["properties"]['dataProvider']['enum'])

    for entry in data["data"]:
        gene_id = entry.get('gene', {}).get('geneId', None)
        if gene_id:
            name, _ = gene_id.split(':', 1)
            assert name in known, "Unknown database: %s" % name

        for global_id in entry['crossReferenceIds']:
            name, _ = global_id.split(':', 1)
            assert name in known, "Xref to unknown db: %s" % name


# Unsure if we should require this, maybe make it an option? I will leave the
# code here for now.
def validate_id_format(data):
    expected = data['metaData']['dataProvider']
    for entry in data['data']:
        primary_id = entry['primaryId']
        db, _ = primary_id.split(':', 1)
        assert db == expected, "Expected %s to start with %s" % (primary_id, expected)

        gene_id = entry.get('gene', {}).get('geneId', None)
        if gene_id:
            gene_db = primary_id.split(':', 1)
            assert gene_db == expected


@click.command()
@click.argument('filename')
@click.option('--schema', default='rnacentral-schema.json',
              help='Filename of the schema to use')
@click.option('--sections', default=SECTIONS,
              help='Directory where schema parts are kept')
def main(filename, schema=None, sections=None):
    with open(filename, 'rb') as raw:
        data = json.load(raw)

    with open(schema, 'rb') as raw:
        schema_data = json.load(raw)

    sections = os.path.abspath(sections)
    base = 'file://%s/' % sections
    file_resolver = js.RefResolver(base_uri=base, referrer=schema)

    js.validate(
        data,
        schema_data,
        format_checker=js.FormatChecker(),
        resolver=file_resolver,
    )
    validate_secondary_structure(data)
    validate_is_known_global_ids(data)


if __name__ == '__main__':
    main()
