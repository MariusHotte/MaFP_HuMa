Dieser Ordner beinhaltet Vorlagen für

a) das python skript zur Bearbeitung eines Versuchs, inklusive custom skripte (zB table.py) im ordner python_custom_scripts
b) latex Protokoll inklusive einzelnen Kapitel files (im Ordner content) sowie header file und einer eigens angelegten mathebefehle.tex für eine vollständigere Bibliothek für mathematische Operatoren (zB partielle Ableitung in schnell)

Ferner sind auch schon Ordnerstrukturen (mit dummys) für messdaten, das fertige Protokoll (final) und ressources (bilder für latex) vorhanden. Die Makefile erlaubt das Kompilieren mit make (mögliche weitere optionen sind 'make once', 'make python', 'make see'). 
Benötigt wird außerdem 
- python 3.4.3 
- lualatex version beta-0.80.0 aus tex live 2015 
python benötigt einige Pakete (numpy, uncertainties, matplotlib..), am besten über Anaconda installieren. Zur Installation der Software siehe 
Pep et al. Toolbox-Workshop 2016.html 
hier im Ordner.
Die Quellen für latex sind in Quellen.bib und programme.bib zu finden und anzupassen. meine-matplotlibrc liefert Einstellungen für schöne Plots mit python mit der Möglichkeit, latex code zu verwenden. Wir haben sehr viele Code-snippets in PythonSkript.py gesammelt, die einem das Leben doch erheblich vereinfachen können. Ein besonderes Augenmerk sei auf die Funktion make_full_table gelegt. Diese erlaubt das einheitliche Erstellen von gut formatierten latex Tabellen aus python heraus. Nun viel Spaß damit ;)