# specset

## Quick start
Create a Python enviroment and install specset

On Windows, we have tested with the official Python distribution Python 3.14.1 at [https://www.python.org/downloads/windows/](https://www.python.org/downloads/windows/). Open a command prompt at a working directory
```
py -m venv venv_specset

venv_specset\Scripts\activate.bat

pip install "specset[full]"
```



On Linux
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

## GUI

There is also an experimental GUI application

```
pip install -e "specset[full]"
```

and start the application by

```
specset-gui
```

One some Linux system, the user has to install a system library

```
sudo apt install libxcb-cursor0
```

or

```
sudo apt install libxcb-xinerama0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-xkb1 libxkbcommon-x11-0
```
