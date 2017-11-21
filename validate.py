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

    # schema_dir = os.path.dirname(os.path.abspath(schema))
    base = 'file://%s/' % sections
    file_resolver = js.RefResolver(base_uri=base, referrer=schema)

    js.validate(
        data,
        schema_data,
        format_checker=js.FormatChecker(),
        resolver=file_resolver,
    )
    validate_secondary_structure(data)


if __name__ == '__main__':
    main()
