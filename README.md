# Processamento de Dados Genômicos e Inferências Filogenéticas
O presente repositório descreve as etapas de processamento de dados genômicos obtidos via *target-enrichment sequencing* com o painel Cactaceae591 (Romeiro-Brito et al. 2022) e a construção de inferências filogenéticas, incluindo a trimagem de dados, montagem de locos, análise de qualidade e construção de árvores gene e espécies. Os programas e comandos necessários devem ser executados no sistema operacional Linux (ou no Subsistema Linux para Windows), através da interface *Bash*.

Base do repositório: Os processos aqui descritos foram adotados no projeto de iniciação científica "Filogenia multilocos de *Arrojadoa* Britton & Rose, cacto endêmico do leste do Brasil", desenvolvido sob orientação do professor Dr. Evandro M. Moraes no Laboratório de Diversidade Genética e Evolução (LaGEvol, UFSCar campus Sorocaba). Os comandos apresentados foram baseados no workshop "Análises Filogenéticas Computacionais", ministrado por M. Kohler, M. Telhe e M. Romeiro-Brito, na mesma instituição. 

## Índice

1. [Filtragem inicial dos dados brutos](#filtragem-inicial-dos-dados-brutos)

2. [Montagem dos locos](#Montagem-dos-locos)

3. [Identificação e remoção de locos possivelmente parálogos](#Identificação-e-remoção-de-locos-possivelmente-parálogos)

4. [Alinhamento de sequências](#Alinhamento-de-sequências)

5. [Polimento das sequências alinhadas](#Polimento-das-sequências-alinhadas)

6. [Inferência filogenética com máxima verossimilhança](#Inferência-filogenética-com-máxima-verossimilhança)

   6.1 [Remoção de regiões hipervariáveis](#Remoção-de-regiões-hipervariáveis)

7. [Inferência filogenética com o método de coalescência](#Inferência-filogenética-com-o-método-de-coalescência)

## Filtragem inicial dos dados brutos
**Objetivo**: Remover sequências de baixa qualidade e adaptadores para garantir uma análise mais precisa dos dados.

Definindo as variáveis referentes às duas fitas do DNA sequenciado:

```
R1=(*R1.fastq)
R2=(*R2.fastq)
```

O programa *fastp* remove sequências de baixa qualidade, muito curtas e adaptadores dos dados. Comandos específicos permitem personalizar valores de parâmetros, mas aqui usaremos os valores padrões fornecidos pelo programa. Para mais informações, consulte: https://github.com/OpenGene/fastp.git.

Filtragem de uma única amostra:

```
fastp -i /caminho/para/a/pasta -I /caminho/para/a/pasta -o /caminho/para/a/pasta/_trimmed.fastq -O /caminho/para/a/pasta/_trimmed.fastq
```

Filtragem de várias amostras: 

```
for ((i=0;i<=${#R1[@]};i++)); do fastp  -i "${R1[i]}" -I "${R2[i]}"  -o "${R1[i]%.fastq}_trimmed.fastq" -O "${R2[i]%.fastq}_trimmed.fastq" -j "${R1[i]}.fastp.json" -h "${R1[i]}.fastp.html" -q 20 --dont_overwrite --failed_out "failed.${R1[i]}"; done
```

Após a execução, os arquivos com os dados filtrados terão o sufixo `_trimmed.fastq`. Para melhor organização, é recomendável movê-los para uma nova pasta. 

## Montagem dos locos
**Objetivo**: Mapear as leituras em locos alvo (assembling), utilizando o *Hybpiper* (https://github.com/mossmatters/HybPiper.git).

Para realizar a montagem dos dados em locos usaremos um arquivo com sequências de referência para cada loco, indicado pela *flag* `-t_dna`.

*Assembling* de uma única amostra:

```
hybpiper assemble -r amostra_trimada_R1.fastq amostra_trimada_R2.fastq -t_dna arquivo_de_referência.fasta --prefix nome_da_amostra --bwa --no_padding_supercontigs --thresh 75
```

*Assembling* de várias amostras:

Crie um arquivo de lista com os nomes das amostras (`namelist.txt`).
```
while read name; do hybpiper assemble -r amostra_trimada_R1.fastq amostra_trimada_R2.fastq -t_dna arquivo_de_referência.fasta --prefix $name --bwa --no_padding_supercontigs --thresh 75; done < namelist.txt
```

O *HybPiper* também possui comandos para gerar estatísticas sobre os locos montados.

Geração de estatísticas do *HybPiper* para os locos montados:

Formato de tabela (`seq_lenghts.tsv`):
```
hybpiper stats -t_dna arquivo_de_referência.fasta gene namelist.txt
```
Formato de mapa de calor:
```
hybpiper recovery_heatmap seq_lengths.tsv
```

Recuperação das leituras de locos em um único arquivo multi-fasta:

```
hybpiper retrieve_sequences dna -t_dna arquivo_de_referência.fasta --sample_names namelist.txt --fasta_dir nome_da_pasta_regioes_alvo
```

Recuperação das leituras de regiões adjacentes aos locos (*supercontigs*):
```
hybpiper retrieve_sequences supercontig -t_dna arquivo_de_referência.fasta --sample_names namelist.txt --fasta_dir nome_da_pasta_supercontigs
```

As estatísticas geradas pelo *HybPiper* podem apontar amostras que não atingiram os valores esperados e que devem ser removidas da base de dados:

Crie um arquivo listando as amostras que deseja remover (`samples_to_delete.txt`).
```
for i in *.fasta; do pxrms -s $i -f /caminho/para/a/pasta/samples_to_delete.txt -o "${i%.FNA}_Reduced.fasta"; done
for i in *.FNA; do pxrms -s $i -f /caminho/para/a/pasta/samples_to_delete.txt -o "${i%.FNA}_Reduced.fasta"; done
```

Remova os avisos adicionados ao nome das suas amostras pelo programa. Utilize os seguintes comandos para os dois formatos de arquivo gerados (fasta e FNA):
```
sed -i 's/\ single_hit//g' *.fasta
sed -i -E 's/\ multi_hit_stitched_contig_comprising_[0-9]+_hits//g' *.fasta
sed -i 's/\.0.*/.0/g' *.fasta
sed -i 's/\.main.*/.main/g' *.fasta
sed -i 's/ .*//' *.fasta
sed -i 's/\./@/g' *.fasta
```

## Identificação e remoção de locos possivelmente parálogos

A detecção de loci parálogos será feita com o *HybPiper* e seguindo os procedimentos de inspeção visual descrito em Frost et al. (2024):

```
hybpiper paralog_retriever namelist.txt -t_dna arquivo_de_referência.fasta
```
Esse comando gera estatísticas sobre os locos parálogos detectados e cria uma pasta chamada `paralogs_all`, contendo as sequências identificadas como parálogas. As sequências desta pasta podem ser alinhadas (com o *MAFFT*) e utilizadas para construção de árvores de genes usando o *IQTree*:
```
nohup sh -c 'for i in *.fasta; do mafft --reorder --auto "$i" > "caminho/para/a/pasta/aligned_$i"; done'  &
iqtree -s aligned_nome_da_amostra_paralogs_all.fasta
```
As árvores geradas devem ser revisadas visualmente para identificar indícios de paralogia, como múltiplas cópias de uma amostra posicionadas de forma distante na árvore. Para orientações adicionais sobre a identificação de parálogos por inspeção visual, acesse Frost et al. (2024): https://doi.org/10.1093/sysbio/syad076. 

![image](https://github.com/user-attachments/assets/d74a76e8-4f6c-4b0c-92f7-2219696fd03a)

Após a identificação, os locos parálogos podem ser removidos da base de dados, movendo-os para uma pasta separada:

Crie uma lista dos locos identificados como parálogos (`paralogs_list.txt`)
```
mkdir locos_paralogos
while read line; do mv $line ./locos_paralogos; done < paralogs_list.txt
```

## Alinhamento de sequências
Para alinhar as sequências de cada loco para todas as amostras, utilizaremos o programa *MAFFT* (https://github.com/GSLBiotech/mafft). É recomendável criar pastas específicas para organizar os alinhamentos. Os arquivos com os alinhamentos terão o prefixo `aligned_`. Nesse repositório, usaremos os parâmetros padrões do programa:

```
nohup sh -c 'for i in *.fasta; do mafft --reorder --auto "$i" > "caminho/para/a/pasta/aligned_$i"; done'  &
```
Os alinhamentos gerados terão o prefixo `aligned_`.

## Polimento das sequências alinhadas
Após o alinhamento, uma segunda etapa de trimagem é necessária para remover posições com alta proporção de *gaps*. Para isso, usaremos o *trimal* (https://github.com/inab/trimal.git). A *flag* `-gt` define o limite de tolerância para *gaps*; por exemplo, `-gt 0.7` remove as "colunas" (sítios) do alinhamento onde a fração de *gaps* é maior ou igual a 30%:

```
mkdir alinhamento_trimado
for i in *.fasta; do trimal -in $i -out ./alinhamento_trimado/"$i"_trimmed.fasta -gt 0.7; done;
```

Para obter estatísticas dos alinhamentos, utilize o *AMAS* (https://github.com/marekborowiec/AMAS.git) para gerar uma planilha com essas informações, o que facilita a análise do impacto de dados faltantes (*gaps*) na base de dados.

```
python3 AMAS.py summary -f fasta -d dna -i *.fasta -o SummaryStats.csv
```

## Inferência filogenética com máxima verossimilhança
Utilizaremos o programa *IQtree* (https://github.com/iqtree/iqtree2.git) para gerar uma árvore de máxima verossimilhança para cada loco:

```
nohup sh -c 'for i in *.fasta; do iqtree2 -nt 4 -s "$i" -st DNA -m MFP -B 10000; done 2>iqtree.err' &
```

Para construção da árvore de máxima verossimilhança para todos os locos, as sequências de cada loco e cada indivíduo precisam ser concatenadas, formando uma supermatriz. A concatenação é feita com o comando:
```
pxcat -s *fasta -p partitions.txt -o minha_supermatriz.fasta
```
A *flag* `-p` gera um arquivo de partições, útil caso seja necessário "desconcatenar" a supermatriz posteriormente.

Uma vez construída a supermatriz, geramos a árvore de espécies:
```
iqtree -nt 4 -s minha_supermatriz.fasta -st DNA -m MFP -B 10000
```
### Remoção de regiões hipervariáveis
Para refinar a supermatriz, é recomendável uma última trimagem com o *spruceup* (https://github.com/marekborowiec/spruceup.git), que remove regiões hipervariáveis (*outliers*). Para executar o programa, é necessário ter a supermatriz, uma árvore de espécies (opcional) e um arquivo de configuração com os parâmetros. Um modelo deste arquivo está anexado ao repositório (`configuration_spruceup.conf`).

```
python -m spruceup configuration-spruceup.conf
```

Esse comando gera uma supermatriz trimada (`0.valor_do_cut_off_minha_supermatriz.fasta`) que pode ser usada para inferir uma nova árvore de espécies pelo método de máxima verossimilhança. A supermatriz trimada também pode ser "desconcatenada" utilizando o arquivo de partições (`partitions.txt`), etapa fundamental para preparar os dados para inferências filogenéticas coalescentes:

```
python3 AMAS.py split -f fasta -d dna -i supermatriz_trimada.fasta -l partitions.txt -u fasta -j
```
Para gerar estatísticas da supermatriz trimada, utilize o AMAS:
```
python3 AMAS.py summary -f fasta -d dna -i *.fasta -o SummaryStats.csv
```

## Inferência filogenética com o método de coalescência
Para construir uma árvore de espécies utilizando o método de coalescência, primeiro é preciso reunir todas as árvores de genes em um único arquivo (all_gene.trees):
```
cat *.treefile > all_gene.trees
```
É recomendável polir este arquivo, reduzindo a influência de nós com baixo suporte. Isso é feito colapsando nós que não atingem um valor de suporte mínimo:

```
nw_ed all_gene.trees 'i & b<=50' o >  all_gene_BS50.tree
```
No comando acima, o valor de suporte mínimo foi definido como 50, mas esse limite pode ser ajustado conforme as características dos dados. Para mais informações sobre esse ajuste, consulte Simmons & Gatesy (2021): https://doi.org/10.1016/j.ympev.2021.107092.

Com as árvores de genes compiladas e polidas, utilizaremos o *ASTRAL III* (https://github.com/smirarab/ASTRAL.git) para inferir a árvore de espécies coalescente:
```
java astral.5.7.8.jar -i all_gene_BS50.tree -o sptree_astral_BS50.tree 2> sptree_astral.log
```

## Referências 

1. FROST, L. A.; BEDOYA, A. M.; LAGOMARSINO, L. P. Artifactual Orthologs and the Need for Diligent Data Exploration in Complex Phylogenomic Datasets: A Museomic Case Study from the Andean Flora. Systematic biology, 3 jan. 2024. DOI: https://doi.org/10.1093/sysbio/syad076

2. ROMEIRO-BRITO, M.; et al. A target Capture Probe Set Useful for Deep - and Shallow - Level Phylogenetic Studies in Cactaceae. Genes, v. 13, n. 707, 2022. DOI: https://doi.org/10.3390/genes13040707.

3. SIMMONS, M. P.; GATESY, J. Collapsing dubiously resolved gene-tree branches in phylogenomic coalescent analyses. Molecular Phylogenetics and Evolution, v. 158, p. 107092, mai. 2021. DOI: https://doi.org/10.1016/j.ympev.2021.107092. 











