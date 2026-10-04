# specset

Create a Python enviroment and install specset

```
python -m venv venv_specset
. venv_specset/bin/activate

pip install specset
```

Follow examples in the notebooks or check the command line usage. For example, the following two command import a MaxQuant search output txt folder to specset and export specset to spectral library

```
specset maxquant_to_specset --additional-column Score=score --additional-column Proteins=proteins --additional-column "Gene Names"=genes txt txt.specset

specset specset_to_speclib --additional-column proteins=Protein.Ids --additional-column genes=Genes --consensus-method score txt.specset txt.specset.lib.tsv
```
