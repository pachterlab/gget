Format of dictionaries containing arguments and expected results:  
```
{
    "test1": {
        "type": "assert_equal",
        "args":{
            "gene": "FUNDC1",
            "sequence": ""
        },
        "expected_result": []
    },
    
    "none_test1": {
        "type": "none",
        "args":{
            "gene": "",
        },
        "expected_result": None
    },
    
    "error_test1": {
        "type": "none",
        "args":{
            "gene": "banana",
        },
        "expected_result": "ValueError"
    },
}
```

Note: Results returned in a DataFrame format need to be converted to a list and NA values should be dropped (df.dropna(axis=1).values.tolist()).  

**Test types:**  
assert_equal -> Use assertListEqual or assertEqual. 
none -> Expected result is None.  
error -> Expected result is specified error.  
