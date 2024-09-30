> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
# gget gpt 💬
Genera texto en lenguaje natural basado en mensaje de entrada. `gget gpt` use la API 'openai.ChatCompletion.create' de [OpenAI](https://openai.com/).
Este módulo, incluido su código, documentación y pruebas unitarias, fue escrito en parte por Chat-GTP3 de OpenAI.

TENGA EN CUENTA:  
Las llamadas a la API de OpenAI solo son 'gratuitas' durante los primeros tres meses después de generar su cuenta de OpenAI (OpenAI proporciona un crédito de $5 que vence).  
Puede definir un límite de facturación mensual estricto (por ejemplo, $1) [aquí](https://platform.openai.com/account/billing/limits).  
Vea sus precios y preguntas frecuentes [aquí](https://openai.com/pricing).  
Obtenga su clave API de OpenAI [aquí](https://platform.openai.com/account/api-keys).  

Regresa: El texto generado (str).  

Antes de usar  `gget gpt` por primera vez, corre `gget setup gpt` / `gget.setup("gpt")` (ver también [`gget setup`](setup.md)).  

**Parámetros posicionales**  
`prompt`  
Mensaje de entrada basado en el cual generar texto (str).  

`api_key`  
Su clave API de OpenAI (str) ([obtenga su clave API](https://platform.openai.com/account/api-keys)).  

**Parámetros optionales**  
`-m` `--model`  
El nombre del algoritmo GPT que se usará para generar el texto (str). Por defecto: "gpt-3.5-turbo".  
See https://platform.openai.com/docs/models/gpt-4 for more information on the available models.  

`-temp` `--temperature`   
Valor entre 0 y 2 que controla el nivel de aleatoriedad y creatividad en el texto generado (float).  
Los valores más altos resultan en un texto más creativo y variado. Por defecto: 1.  

`-tp` `--top_p`   
Controla la diversidad del texto generado como alternativa al muestreo con `--temperature` (float).  
Los valores más altos resultan en un texto más diverso e inesperado. Por defecto: 1.  
Tenga en cuenta que OpenAI recomienda modificar `--top_p` o el parámetro `--temperature`, pero no ambas.  

`-s` `--stop`   
Una secuencia de tokens para marcar el final del texto generado (str). Por defecto: None.  

`-mt` `--max_tokens`   
Controla la longitud máxima del texto generado, en tokens (int). Por defecto: 200.  

`-pp` `--presence_penalty`   
Número entre -2.0 y 2.0. Los valores más altos aumentan la probabilidad de que el modelo hable sobre temas nuevos (float). Por defecto: 0.  

`-fp` `--frequency_penalty`   
Número entre -2.0 y 2.0. Los valores más altos reducen la probabilidad de que el modelo repita la misma línea palabra por palabra (float). Por defecto: 0.  

`-lb` `--logit_bias`   
Un diccionario que especifica un sesgo hacia ciertos tokens en el texto generado (dict). Por defecto: None.  

`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados, p. ej. ruta/al/directorio/resultados.txt. Por defecto: salida estándar (STDOUT).  
  
### Por ejemplo
```bash
gget gpt "Cómo estás hoy GPT?" su_clave_api
```
```python
# Python
print(gget.gpt("Cómo estás hoy GPT?", "su_clave_api"))
```

<br>

<img width="725" alt="Screen Shot 2023-03-18 at 3 42 32 PM" src="https://user-images.githubusercontent.com/56094636/226143902-6fa2d0c7-7eea-4382-b1d2-df6c3f0d5fd5.png">
