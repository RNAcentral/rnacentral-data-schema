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

import pytest

import validate

EXAMPLE_DIR = 'examples'


@pytest.mark.parametrize(
    'filename',
    [os.path.join(EXAMPLE_DIR, f) for f in os.listdir(EXAMPLE_DIR)]
)
def test_can_validate_all_example_files(filename):
    with open(filename, 'r') as raw:
        data = json.load(raw)

    print(os.path.abspath(os.curdir))
    print(os.listdir('sections'))
    validate.validate(data, validate.SCHEMA_NAME)
