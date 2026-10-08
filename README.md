# specset

## Installation
Create a Python enviroment and install specset

On Windows, we have tested with the official Python distribution Python 3.14.1 at [https://www.python.org/downloads/windows/](https://www.python.org/downloads/windows/). Open a command prompt at a working directory
```
py -m venv venv_specset

venv_specset\Scripts\activate.bat

pip install "specset[full]"
```
The ```full``` installation includes a GUI application. 


On Linux, the installation is similar
```
python -m venv venv_specset

. venv_specset/bin/activate

pip install "specset[full]"
```

For the GUI application, one some Linux systems, the user has to install extra system libraries

```
sudo apt install libxcb-cursor0
```

or

```
sudo apt install libxcb-xinerama0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-xkb1 libxkbcommon-x11-0
```


For experiments with AlphaPeptDeep, installation of the AlphaPeptDeep ecosystem is needed.
```
pip install "alphabase[full-stable]"
pip install "peptdeep[stable,gui-stable]"
```

## Usage

Follow examples in the notebooks or check the command line usage. For example, the following two commands import a MaxQuant search output txt folder into a specset and export the specset to a spectral library

```
specset maxquant_to_specset --additional-column Score=score --additional-column Proteins=proteins --additional-column "Gene Names"=genes txt txt.specset

specset specset_to_speclib --additional-column proteins=Protein.Ids --additional-column genes=Genes --consensus-method score txt.specset txt.specset.lib.tsv
```

## GUI application

The application can be started by executing the command
```
specset-gui
```
The application features three tabs: Viewer, Compare, and Tools. The Viewer tab allows users to browse a single spec set. The Compare tab highlights the differences between two spec sets. The Tools tab directly maps to the available shell commands.

