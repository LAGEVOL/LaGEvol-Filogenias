# Processamento de Dados Genômicos e Inferências Filogenéticas
Repositório sobre o processamento de dados genômicos e inferências filogenéticas. Os dados utilizados como modelo foram sequenciados utilizando o painel de captura Cactaceae591. 

## Trimagem inicial dos dados brutos
Definindo as variáveis referentes às duas fitas do DNA sequenciado.

```
R1=(*R1.fastq)
R2=(*R2.fastq)
```

Para esta trimagem inicial dos dados brutos, utilizaremos o *fastp*. Este programa remove sequências de baixa qualidade, sequências muito curtas e adaptadores dos nossos dados. Estes parâmetros podem ser modificados com comandos específicos, mas neste caso utilizaremos os valores padrões fornecidos pelo programa. Para mais informações, acesse: https://github.com/OpenGene/fastp.git.

Utilizando o *fastp* para uma única amostra:

```
fastp -i /caminho/para/a/pasta -I /caminho/para/a/pasta -o /caminho/para/a/pasta/_trimmed.fastq -O /caminho/para/a/pasta/_trimmed.fastq
```

Utilizando o *fastp* para várias amostras: 

```
for ((i=0;i<=${#R1[@]};i++)); do fastp  -i "${R1[i]}" -I "${R2[i]}"  -o "${R1[i]%.fastq}_trimmed.fastq" -O "${R2[i]%.fastq}_trimmed.fastq" -j "${R1[i]}.fastp.json" -h "${R1[i]}.fastp.html" -q 20 --dont_overwrite --failed_out "failed.${R1[i]}"; done
```

Os arquivos com os dados trimados terão o sufixo `_trimmed.fastq`. Para melhor organização, é recomendável que estes arquivos sejam movidos para uma nova pasta. 


