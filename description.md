# A lightweight toolkit for manipulation of peptide spectrum matches in mass spectrometry-based proteomics

## Installation
Create a Python enviroment and install specset using pip

### Windows

On Windows, we have tested with the official Python distribution Python 3.12.9 at [https://www.python.org/downloads/windows/](https://www.python.org/downloads/windows/). Open a command prompt from a working directory
```
python -m venv venv_specset

venv_specset\Scripts\activate.bat

pip install specset[full]
```
The ```full``` installation includes a GUI application. 

### Linux
On Linux, the installation is as follows.
```
python -m venv venv_specset

. venv_specset/bin/activate

pip install specset[full]
```

For the GUI application, one some Linux systems, the user has to install extra system libraries

```
sudo apt install libxcb-cursor0
```

or

```
sudo apt install libxcb-xinerama0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-xkb1 libxkbcommon-x11-0
```


## Usage

Check the command line usage.

```
specset --help
```

The GUI application can be started by executing the command

```
specset-gui
```
