import unittest
from unittest.mock import patch
from gget.gget_gpt import gpt


class TestGpt(unittest.TestCase):
    @patch("openai.Completion.create")
    def test_gpt(self, mock_create):
        # Mock the response from the API
        mock_response = {"choices": [{"text": "This is a generated response."}]}
        mock_create.return_value = mock_response

        # Call the function with a prompt and test API key
        prompt = "Test prompt"
        api_key = "test_api_key"
        output = gpt(prompt, api_key)

        # Check that the mock API was called with the correct arguments
        mock_create.assert_called_once_with(
            engine="davinci",
            prompt=prompt,
            max_tokens=1024,
            stop=None,
            temperature=0.5,
            api_key=api_key,
        )

        # Check that the output matches the mock response
        expected_output = "This is a generated response."
        self.assertEqual(output, expected_output)
