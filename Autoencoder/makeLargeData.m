data0 = load('AML_genomes_mutational_catalog_96_subs.mat');
data1 = load('Breast_genomes_mutational_catalog_96_subs.mat');
data2 = load('CLL_genomes_mutational_catalog_96_subs.mat');
data3 = load('Liver_genomes_mutational_catalog_96_subs.mat');
data4 = load('Lung Adeno_genomes_mutational_catalog_96_subs.mat');
data5 = load('Lymphoma B-cell_genomes_mutational_catalog_96_subs.mat');
data6 = load('Medulloblastoma_genomes_mutational_catalog_96_subs.mat');
data7 = load('Pancreas_genomes_mutational_catalog_96_subs.mat');
data8 = load('Pilocytic Astrocytoma_genomes_mutational_catalog_96_subs.mat');

originalGenomes = cat(2, data0.originalGenomes, data1.originalGenomes, ...
    data2.originalGenomes, data3.originalGenomes, ...
    data4.originalGenomes, data5.originalGenomes, ...
    data6.originalGenomes, data7.originalGenomes, ...
    data8.originalGenomes);

destination = strcat('C:\Users\....mat');
save(destination, 'originalGenomes', 'subtypes', 'types');
