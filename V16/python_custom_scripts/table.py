import itertools
import codecs
import uncertainties
import fileinput
import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds,
)
from uncertainties import ufloat

def search_replace_within_file(filenameToSearch, textToSearch, textToReplace):
    """Function to search files for s sequency and replace this one by a given value
    Args:
            filenameToSearch:   (str) Name of the file
            textToSearch:       (str) Which sequency are you looking for?
            textToReplace:      (str) Replace hits by this sequency
    """
    pass
    with fileinput.FileInput(filenameToSearch, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(textToSearch, textToReplace), end='')

def make_table(columns, figures=None):
    if figures is None:
        figures = [None] * len(columns)

    cols = []
    for column, figure in zip(columns, figures):
        if (type(column) == str):
            column = [column]
        if (type(column) == list):
            col = zip(*zip(column))     # hard to find this kind of code... this will unzip the list column, ref: https://docs.python.org/3/library/functions.html#zip
            cols.extend(col)
        elif np.any(stds(column)):
            if figure is None:
                figure = ''
            col = list(zip(*['{0:.{1:}uf}'.format(x, figure).split('+/-') for x in column]))
        else:
            try:
                test_iterator = iter(column)    # if only one entry is given, then this will throw a type error exception - handled below
                col = list(zip(*[['{0:.{1:}f}'.format(x, figure)] for x in noms(column)]))
            except TypeError:
                col = list(zip(*[['{0:.{1:}f}'.format(column, figure)]]))
        cols.extend(col)

    max_lens = [max(len(s) for s in col) for col in cols]
    cols = [['{0:<{1:}}'.format(s, ml) for s in col] for col, ml in zip(cols, max_lens)]

    rows = list(itertools.zip_longest(*cols))

    return (r' \\' + '\n').join([' & '.join(s for s in row if s is not None) for row in rows]) + r' \\'

def make_composed_table(tables):
    assert isinstance(tables, list), "You need to give a list of filenames to make_composed_table!"
    Output = ''
    for filename in tables:
        with open(filename, 'r') as f:
            Output += f.read()
    return Output

def make_SI(num, unit, exp='', figures=None):
    y = ufloat(0.0, 0) #siunitx mag kein 0 +- 0, deshalb hier der workaround
    if num == y:
        return "(0 \pm 0) ~ \si{" + unit + "}"
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

def write(filename, content):
    f = codecs.open(filename, "w", "utf-8")
    test_type = ufloat(2,3)
    if type(content) == type(test_type):
        content = "\num{" + str(x.n) + " +- " + str(x.s) + "}"
        f.write(content)
        if not content.endswith('\n'):
            f.write('\n')
        f.close()
    else:
        f.write(content)
        if not content.endswith('\n'):
            f.write('\n')
        f.close()

    #
    # with open(filename, 'w') as f:
    #     f.write(content)
    #     if not content.endswith('\n'):
    #         f.write('\n')


def make_full_table(caption,label,source_table, stacking=np.array([]), units=None, replaceNaN = False, replaceNaNby = '-'):
    """This function will create a complete tex-table for your usage. It is formated automatically so you can't choose any format designs.
    Args:
            caption:        (str) The caption which will be shown in the .pdf generated by TeX
            label:          (str) The label which you need for \ref within TeX
            source_table    (str) A valid filename for the source file. Generate it using make_table. You can also use columns of strings.
            stacking        (list) A list with the numbers of columns which will be connected as a consequence of mean-value +- error. Default: None
            units           (list) A list for every resulting column name. Default: None
            replaceNaN      (bool) Set it to true if you want to replace NaN-values in the source file with the parameter replaceNaNby. Default: False
            replaceNaNby    (str) Choose what you want the NaN's from your sourcefile replaced by. Default: '-'
    """
    pass
    # Zunächstmal behandeln wir die Quelldatei falls das gewünscht ist, damit NaN's mit '-' angezeigt werden oder was auch immer man eingegeben hat
    if (replaceNaN):
        search_replace_within_file(source_table, 'nan', replaceNaNby)
    # Zuerst definieren wir eine Funktion, die erkennt, ob es sich bei einem char um einen integer handelt:
    def RepresentsInt(s):
        try:
            int(s)
            return True
        except ValueError:
            return False


    # Vorgeplänkel
    Output = """\\begin{table}
    \\centering
    \\caption{""" + caption + """}
    \\label{""" + label + """}
    \\sisetup{parse-numbers=false}
    \\begin{tabular}{\n"""

    # Kerngeschäft : source_table einlesen und verarbeiten, dh. Vor und Nachkommastellen rausfinden
    counter_columns = 0
    counter_lines = 0
    with open(source_table, 'r') as f:
        Text = f.read()
        for buchstabe in Text:
            if (buchstabe == '&'):
                counter_columns += 1
            elif (buchstabe == '\\'):
                counter_lines += 1

    NumberOfLines = int(counter_lines/2)
    NumberOfColumns = int(counter_columns/counter_lines*2+1)
    counter_digits_preDot = np.zeros((NumberOfLines, NumberOfColumns), dtype=np.int)
    counter_digits_postDot = np.zeros((NumberOfLines, NumberOfColumns), dtype=np.int)
    dot_reached = False
    counter_columns = 0
    counter_lines = 0
    is_last_column_a_string = True
    remember_columns_with_strings = []
    with open(source_table, 'r') as f:
        Text = f.read()
    # 'Vor und Nachkommastellen rausfinden' beginnt hier
        for buchstabe in Text:
            if (buchstabe == '&'):
                if (is_last_column_a_string & (counter_lines==0)):
                    remember_columns_with_strings.append(counter_columns)
                counter_columns += 1
                dot_reached = False
                is_last_column_a_string = True  # wir gehen davon aus, dass ein string drin steht, und ändern diese Meinung wenn wir ne Zahl finden (s.u.),
                                                # das ist also auch wahr falls NICHTS in der spalte stand
            elif (buchstabe == '.'):
                dot_reached = True
            elif (buchstabe == '\\'):
                if (is_last_column_a_string & (counter_lines==0)):
                    remember_columns_with_strings.append(counter_columns)
                counter_lines += 1
                counter_columns = counter_columns % (NumberOfColumns-1)
                dot_reached = False
                is_last_column_a_string = True  # wir gehen davon aus, dass ein string drin steht, und ändern diese Meinung wenn wir ne Zahl finden (s.u.)
            elif (RepresentsInt(buchstabe)):
                is_last_column_a_string = False     # sobald wir auch nur irgendeine Zahl finden ist dieser Wert halt falsch
                if (counter_lines/2 <= (NumberOfLines-1)):
                    if dot_reached == False:
                        counter_digits_preDot[int(counter_lines/2)][counter_columns] += 1
                    else:
                        counter_digits_postDot[int(counter_lines/2)][counter_columns] += 1
    # jetzt ermittle maximale Anzahl an Stellen und speichere sie in MaxDigitsPreDot und MaxDigitsPostDot
    MaxDigitsPreDot = []
    counter_digits_preDot_np = np.array(counter_digits_preDot)
    for x in counter_digits_preDot_np.T:
        MaxDigitsPreDot.append(max(x))
    MaxDigitsPostDot = []
    counter_digits_postDot_np = np.array(counter_digits_postDot)
    for x in counter_digits_postDot_np.T:
        MaxDigitsPostDot.append(max(x))
    # --------------------Ende der Stellensuche

    # Die Liste stacking in ein angepasstes Array umwandeln mit den tatsächlich betroffenen Spalten
    stacking_list = np.array(stacking)
    i = 0
    for x in stacking_list:
        stacking_list[i] += i
        i += 1

    # Schreiben der Tabellenformatierung
    if np.size(stacking) == 0:
        i = 0.0
        for digits_preDot, digits_postDot in zip(MaxDigitsPreDot, MaxDigitsPostDot):
            if i in remember_columns_with_strings:
                Output += '\tl\n'       # left aligned strings
            else:
                Output += '\tS[table-format=' + str(digits_preDot) + '.' + str(digits_postDot) +']\n'
            i += 1
    else:   # es wurden fehlerbehaftete Werte übergeben, daher muss +- zwischen die entsprechenden Spalten
        i = 0.0
        for digits_preDot, digits_postDot in zip(MaxDigitsPreDot, MaxDigitsPostDot):
            if i in remember_columns_with_strings:
                Output += '\tl\n'       # left aligned strings
            elif i in stacking_list:
                Output += '\tS[table-format=' + str(digits_preDot) + '.' + str(digits_postDot) +']\n'
                Output += '\t@{${}\\pm{}$}\n'
            elif i-1 in stacking_list:
                Output += '\tS[table-format=' + str(digits_preDot) + '.' + str(digits_postDot) +', table-number-alignment = left]\n'      # wir wollen hier linksbündige Zahlen
            else:
                Output += '\tS[table-format=' + str(digits_preDot) + '.' + str(digits_postDot) +']\n'
            i += 1

    # Zwischengeplänkel
    Output += '\t}\n\t\\toprule\n\t'

    # Einheitenzeile
    i=0
    stacking_list = np.array(stacking)
    for Spaltenkopf in units:
        if i in stacking_list:
            Output += '\\multicolumn{2}{c}'
        Output += '{' + str(Spaltenkopf) + '}\t\t'
        i += 1
        if i == np.size(units):
            Output += '\\\\ \n\t'
        elif i % 2 == 0:
            Output += '& \n\t'
        else:
            Output += '& '

    # Schlussgeplänkel
    Output += """\\midrule
    \\input{""" + source_table + """}
    \\bottomrule
    \\end{tabular}
    \\end{table}"""
    return Output
