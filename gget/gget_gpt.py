import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

import openai


def gpt(
    prompt,
    api_key,
    engine="davinci",
    max_tokens=100,
    stop=None,
    temperature=0.5,
    output=None,
):
    """
    Generates text based on a given prompt using the OpenAI API.

    Args:
        prompt (str): The input prompt for the GPT-3 model to generate text from.
        api_key (str): Your OpenAI API key (see: https://platform.openai.com/account/api-keys).
        engine (str, optional): The name of the GPT-3 engine to use (defaults to "davinci").
                                You can choose from the following engines:
                                "davinci": Most capable, creative, and accurate but slower.
                                "curie": A good balance between performance and response time.
                                "babbage": Faster and lower-cost with decent performance.
                                "ada": Fastest and least expensive, suitable for simple tasks.
        max_tokens (int, optional): The maximum number of tokens (words or subwords) in the generated text (defaults to 100).
        stop (str, optional): A sequence of tokens that should indicate the end of the generated text (defaults to None).
        temperature (float, optional): Controls the "creativity" of the generated text (defaults to 0.5).

    Returns:
        str: The generated text based on the input prompt.

    Raises:
        openai.Error: If there is an error with the OpenAI API request.

    Example:
        >>> prompt = "What is the meaning of life?"
        >>> api_key = "YOUR_API_KEY_HERE"
        >>> import gget
        >>> response = gget.gpt(prompt, api_key, engine="curie", max_tokens=64, stop="\n", temperature=0.7)
        >>> print(response)
        "The meaning of life is subjective and varies from person to person. Some people find meaning in their work or career, while others find it in their relationships or hobbies."

    NOTE: OpenAI API calls are only 'free' for the first three months after generating your OpenAI Account
    (they give you a $5 credit that expires).
    You can define a hard billing limit (e.g. $1) here: https://platform.openai.com/account/billing/limits
    See their pricing and FAQ here: https://openai.com/pricing
    You can get your API key here: https://platform.openai.com/account/api-keys

    This module, including its source code, documentation and unittests, were written (almost) entirely by OpenAI's Chat-GTP3.
    """

    engines = ["davinci", "curie", "babbage", "ada"]
    if engine not in engines:
        raise ValueError(
            f"Engine specified is {engine}. Expected one of: {', '.join(engines)}"
        )

    openai.api_key = api_key

    response = openai.Completion.create(
        engine=engine,
        prompt=prompt,
        max_tokens=max_tokens,
        n=1,
        stop=stop,
        temperature=temperature,
    )

    logging.info(
        f"Total tokens used for API call to engine '{engine}': {response['usage']['total_tokens']}"
    )

    texts = response["choices"][0]["text"].strip()

    if output:
        with open(output, "w") as f:
            f.write(texts)

    return texts
