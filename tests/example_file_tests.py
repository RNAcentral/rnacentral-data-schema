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

import copy
import json
import itertools as it
from pathlib import Path

import pytest

import validate


def examples():
    data = []
    for example in Path("examples/").glob("*.json"):
        with example.open("r") as raw:
            data.append(json.load(raw))
    return data


def filenames():
    data = []
    for example in Path("examples/").glob("*.json"):
        data.append(str(example))
    return data


@pytest.mark.parametrize("data", examples())
def test_can_validate_all_example_files(data):
    validate.validate(data, validate.SCHEMA_NAME, validate.SECTIONS)


@pytest.mark.parametrize(
    "field,data", list(it.product(["primaryId", "sequence", "url"], examples()))
)
def test_fails_when_missing_key_properties(field, data):
    del data["data"][0][field]
    with pytest.raises(Exception):
        validate.validate(data, validate.SCHEMA_NAME, validate.SECTIONS)


@pytest.mark.parametrize("filename", filenames())
def test_requires_one_of_taxon_id_or_inferred(filename):
    with open(filename, "r") as raw:
        data = json.load(raw)
    validate.validate(data, validate.SCHEMA_NAME, validate.SECTIONS)
    fields = ["taxonId", "inferredPhylogeny"]
    for index, ncrna in enumerate(data["data"]):
        current = copy.deepcopy(data)
        for field in fields:
            if field in ncrna:
                del current["data"][index][field]
                break
        else:
            raise ValueError("Must update test to handle no taxon/inferredPhylogeny")

        with pytest.raises(Exception):
            validate.validate(current, validate.SCHEMA_NAME, validate.SECTIONS)


@pytest.mark.parametrize("filename", filenames())
def test_requires_one_of_so_term_or_rfam_accession(filename):
    with open(filename, "r") as raw:
        data = json.load(raw)
    validate.validate(data, validate.SCHEMA_NAME, validate.SECTIONS)
    fields = ["soTermId", "rfamAccession"]
    for index, ncrna in enumerate(data["data"]):
        current = copy.deepcopy(data)
        for field in fields:
            if field in ncrna:
                del current["data"][index][field]
                break
        else:
            raise ValueError("Must update test to handle no taxon/inferredPhylogeny")

        with pytest.raises(Exception):
            validate.validate(current, validate.SCHEMA_NAME, validate.SECTIONS)
