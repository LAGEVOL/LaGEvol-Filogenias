# Processamento de Dados Genômicos e Inferências Filogenéticas
Repositório sobre o processamento de dados genômicos e inferências filogenéticas. Os dados utilizados como modelo foram sequenciados utilizando o painel de captura Cactaceae591. 

## Trimagem inicial dos dados brutos
Definindo as variáveis referentes às duas fitas do DNA sequenciado:

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

## Montagem dos dados trimados em locos
Para montar os dados trimados em seus respectivos locos (processo chamado de *assembling*), utilizaremos o programa *HybPiper* (https://github.com/mossmatters/HybPiper.git). Precisaremos de um arquivo de referência para a montagem dos locos, indicado nos comandos pela *flag* `-t_dna`. 

Utilizando o *HybPiper* em uma única amostra:

```
hybpiper assemble -r amostra_trimada_R1.fastq amostra_trimada_R2.fastq -t_dna arquivo_de_referência.fasta --prefix nome_da_amostra --bwa --no_padding_supercontigs --thresh 75
```

Utilizando o *HybPiper* em várias amostras:

Primeiramente, crie um arquivo listando as suas amostras (ex: namelist.txt).
```
while read name; do hybpiper assemble -r amostra_trimada_R1.fastq amostra_trimada_R2.fastq -t_dna arquivo_de_referência.fasta --prefix $name --bwa --no_padding_supercontigs --thresh 75; done < namelist.txt
```

O *HybPiper* também possui comandos para gerar estatísticas sobre os locos montados.

Gerando estatísticas no formato de tabela (seq_lengths.tsv):
```
hybpiper stats -t_dna arquivo_de_referência.fasta gene namelist.txt
```
Gerando estatísticas no formato de mapa de calor:
```
hybpiper recovery_heatmap seq_lengths.tsv
```

Os *reads* recuperados para cada região alvo devem ser sintetizados em um único arquivo multi-fasta, com o comando:

```
hybpiper retrieve_sequences dna -t_dna arquivo_de_referência.fasta --sample_names namelist.txt --fasta_dir nome_da_pasta_regioes_alvo
```

Utilizando o *HybPiper*, também é possível recuperar regiões adjacentes às regiões alvo do painel de captura utilizado no sequenciamento. A combinação de regiões alvo do painel + sequências adjacentes é chamada de *supercontig*, e pode ser de grande valor a depender dos dados e análises trabalhadas.
```
hybpiper retrieve_sequences supercontig -t_dna arquivo_de_referência.fasta --sample_names namelist.txt --fasta_dir nome_da_pasta_supercontigs
```

Analisando as estatísticas geradas pelo *HybPiper*, você pode encontrar amostras que não atingiram os parâmetros esperados. Neste caso, você pode removê-las da sua base de dados:

Crie um arquivo listando as amostras que deseja remover (ex: samples_to_delete.txt).
```
for i in *.fasta; do pxrms -s $i -f /caminho/para/a/pasta/samples_to_delete.txt -o "${i%.FNA}_Reduced.fasta"; done
for i in *.FNA; do pxrms -s $i -f /caminho/para/a/pasta/samples_to_delete.txt -o "${i%.FNA}_Reduced.fasta"; done
```

Por fim, é importante que você remova avisos que o programa pode ter adicionado ao nome das suas amostras. Para isso, utilize os seguintes comandos para os dois formatos de arquivo gerados (fasta e FNA):
```
sed -i 's/\ single_hit//g' *.fasta
sed -i -E 's/\ multi_hit_stitched_contig_comprising_[0-9]+_hits//g' *.fasta
sed -i 's/\.0.*/.0/g' *.fasta
sed -i 's/\.main.*/.main/g' *.fasta
sed -i 's/ .*//' *.fasta
sed -i 's/\./@/g' *.fasta
```

## Identificação e remoção de locos possivelmente parálogos

## Alinhamento de sequências

## Polimento das sequências alinhadas








