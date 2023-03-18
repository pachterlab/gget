> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget gpt 💬
Generates text based on a given prompt using the [OpenAI](https://openai.com/) API.  
This module, including its source code, documentation and unit tests, were written (almost) entirely by OpenAI's Chat-GTP3.  

NOTE:  
OpenAI API calls are only 'free' for the first three months after generating your OpenAI Account (OpenAI provides a $5 credit that expires).  
You can define a hard monthly billing limit (e.g. $1) [here](https://platform.openai.com/account/billing/limits).  
See their pricing and FAQ [here](https://openai.com/pricing).  
Get your OpenAI API key [here](https://platform.openai.com/account/api-keys).  

Returns: Predicted API response (text).  

**Positional argument**  
`prompt`  
The input prompt for the GPT-3 model to generate text from (str).  

`api_key`  
Your OpenAI API key (str) ([get your API key](https://platform.openai.com/account/api-keys)).  

**Optional arguments**  
`-e` `--engine`  
The name of the GPT-3 engine to use (defaults to 'davinci').  
You can choose from the following engines:  
    "davinci": Most capable, creative, and accurate but slower.  
    "curie": A good balance between performance and response time.  
    "babbage": Faster and lower-cost with decent performance.  
    "ada": Fastest and least expensive, suitable for simple tasks.  

`-m` `--max_tokens`   
The maximum number of tokens (words or subwords) in the generated text (defaults to 100).  

`-s` `--stop`   
A sequence of tokens that should indicate the end of the generated text (defaults to None).  

`-temp` `--temperature`   
Controls the 'creativity' of the generated text (defaults to 0.5).  
A higher value increases creativity, while a lower value decreases creativity.  

`-o` `--output`   
The file name to save the generated text to as a text file (defaults to printing the output to the console).  
  
  
### Example
```bash
gget gpt "hi there gpt" your_api_token
```
```python
# Python
gget.gpt("hi there gpt", "your_api_token")
```

<br>

<img width="725" alt="Screen Shot 2023-03-18 at 3 42 32 PM" src="https://user-images.githubusercontent.com/56094636/226143902-6fa2d0c7-7eea-4382-b1d2-df6c3f0d5fd5.png">