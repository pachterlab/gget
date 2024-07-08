import unittest
import json
from typing import Literal
from gget.gget_bgee import bgee_expression, bgee_orthologs
from .from_json import from_json

# Load dictionary containing arguments and expected results
with open("./tests/fixtures/test_bgee.json") as json_file:
    bgee_dict = json.load(json_file)


# noinspection PyShadowingBuiltins
def fun(type: Literal["expression", "orthologs"], *args, **kwargs):
    if type == "expression":
        return bgee_expression(*args, **kwargs)
    elif type == "orthologs":
        return bgee_orthologs(*args, **kwargs)
    else:
        raise ValueError(f"Unknown type: {type}")


class TestBgee(unittest.TestCase, metaclass=from_json(bgee_dict, fun)):
    pass  # all tests are loaded from json
