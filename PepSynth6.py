import sys
from os import path, listdir
import requests

from PyQt6.QtGui import QIcon, QDoubleValidator,QRegularExpressionValidator,QPixmap,QDragEnterEvent
from PyQt6.QtCore import Qt, QAbstractTableModel, QRegularExpression
from PyQt6.QtWidgets import QWidget, QApplication, QStackedWidget, QListWidget, QAbstractItemView, QFileDialog, QHBoxLayout, QComboBox,QScrollArea, QFrame, QGridLayout, QGroupBox
from PyQt6.QtWidgets import QFormLayout, QSpinBox,QPushButton, QLabel,QLineEdit,QMessageBox, QTableView, QVBoxLayout,QCheckBox, QSlider, QSlider, QTextEdit
from qt_material import apply_stylesheet
from PyQt6.QtWebEngineWidgets import QWebEngineView

import pandas as pd
import numpy as np
from matplotlib.figure import Figure
import datetime
import traceback
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from seaborn import scatterplot,move_legend
import re

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=10, height=4, dpi=100, nrows=1, ncols=2):
        self.fig, self.ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)

class TableModel(QAbstractTableModel):
    def __init__(self, data):
        super().__init__()
        self._data = data.copy()

    def data(self, index, role):
        if role == Qt.ItemDataRole.DisplayRole:
            value = self._data.iloc[index.row(), index.column()]
            return str(value)

    def rowCount(self, index):
        return self._data.shape[0]

    def columnCount(self, index):
        return self._data.shape[1]

    def headerData(self, section, orientation, role):
        if role == Qt.ItemDataRole.DisplayRole:
            if orientation == Qt.Orientation.Horizontal:
                return str(self._data.columns[section])

            if orientation == Qt.Orientation.Vertical:
                return str(self._data.index[section])

class PepSyn(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.miniCardNr = 1
        self.pepToStoreVal = 20
        self.stack1Layout()
        self.stack2Layout()
        self.stack3Layout()
        self.stack4Layout()
        self.pepCtl = ["QPAVRNERA", "QPDGGQPAV", "QPAVRNERAA", "PDGGQPAVR", "SFFSYGEIP", "SDSSFFSYG", "SSFFSYGEI",
                        "FSYGEIPFG", "ACENIAECRTSPRCA", "ACEKIAECRTSPRCA", "NVARQPARE", "RRNAAVAEQP", "IPYSGFESF",
                        "FIFYSEGSS", "RTACCINAESECPAR"]
        self.infoText = ''

    def initUI(self):
        self.resize(1320, 800)
        self.setWindowTitle('Demo PepDesigner')
        overLayout = QHBoxLayout()
        vLayout = QVBoxLayout()

        self.stack = QStackedWidget()
        self.stack1 = QWidget()
        self.stack2 = QWidget()
        self.stack3 = QWidget()
        self.stack4 = QWidget()
        self.stack.addWidget(self.stack1)
        self.stack.addWidget(self.stack2)
        self.stack.addWidget(self.stack3)
        self.stack.addWidget(self.stack4)
        self.list = QListWidget()
        self.list.addItems(['Data Input','Peptide Generation', 'Minicard Generation', 'Array Generation'])
        self.list.setCurrentRow(0)
        self.list.currentRowChanged.connect(self.display)

        vLayout.addWidget(self.list)
        overLayout.addLayout(vLayout,2)
        overLayout.addWidget(self.stack,7)

        self.setLayout(overLayout)


    def display(self, index):
        self.stack.setCurrentIndex(index)

    def stack1Layout(self):
        layoutStack1 = QVBoxLayout()
        zeroRowLayout = QGridLayout()
        zeroRow_Group = QGroupBox(self)
        zeroRow_Group.setTitle("&Parameters")
        zeroRow_Group.setLayout(zeroRowLayout)

        thirdRowLayout = QHBoxLayout()
        pdbLayout = QHBoxLayout()
        pdb_Group = QGroupBox(self)
        pdb_Group.setTitle("&Conformational structure")
        pdb_Group.setLayout(pdbLayout)

        lenLayout = QGridLayout()
        lenGroup = QGroupBox(self)
        lenGroup.setTitle("&Length")
        lenGroup.setLayout(lenLayout)

        len1Layout = QHBoxLayout()
        len2Layout = QHBoxLayout()
        len3Layout = QHBoxLayout()
        len4Layout = QHBoxLayout()
        len5Layout = QHBoxLayout()
        len6Layout = QHBoxLayout()
        btLayout = QHBoxLayout()

        # User info
        self.userLabel = QLabel("User:")
        self.userCombo = QComboBox(self)
        self.userCombo.setPlaceholderText("select user")
        self.userCombo.addItems(['Fahui Liu', 'Fahui Liu', 'Fahui Liu',
                                 'Fahui Liu', 'Fahui Liu', 'Fahui Liu'])
        self.userCombo.setFixedWidth(160)

        self.projNameLabel = QLabel(" Project name:")
        self.projNameText = QLineEdit()
        self.projNameText.setFixedWidth(160)

        self.synNumLabel = QLabel("Synthesis number:")
        self.synNumText = QLineEdit()
        self.synNumText.setFixedWidth(160)

        self.synTypeLabel = QLabel("Synthesis type:")
        self.synTypeText = QComboBox(self)
        self.synTypeText.setPlaceholderText("select synth type")
        self.synTypeText.addItems(['LIN','LOOP','HEL','BET','MAT.LIN','MAT.LOOP'])
        self.synTypeText.setFixedWidth(160)

        # target info
        self.targetLabel = QLabel("Target:")
        self.targetText = QLineEdit()
        self.targetText.setFixedWidth(160)

        self.pdbLabel = QLabel("TargetPDB:")
        self.pdbText = QLineEdit()
        self.pdbText.setPlaceholderText("4-letter ID from PDB database")
        self.pdbText.editingFinished.connect(self.pdbTextChanged)
        self.pdbText.setFixedWidth(160)
        self.nSeqLabel = QLabel("n Lengths:")
        self.nSeqText = QSpinBox()
        self.nSeqText.setFixedWidth(160)
        self.nSeqText.setMinimum(1)
        self.nSeqText.setMaximum(6)
        self.nSeqText.valueChanged.connect(self.nLengthChanged)
        self.seqNrLabel = QLabel("Sequence code:")
        self.seqNrText = QSpinBox()
        self.seqNrText.setFixedWidth(160)
        self.seqNrText.setValue(1)
        self.seqNrText.setMinimum(1)
        self.seqNrText.setMaximum(10)

        self.seqLabel = QLabel("Sequence:")
        self.seqLenText = QLabel('Len')
        self.seqText = QLineEdit()
        self.seqText.setPlaceholderText("Paste or type your protein sequence")
        self.seqText.setMaxLength(700)
        self.seqText.editingFinished.connect(self.seqPasted)
        self.seqLen = 20

        floatValidator = QRegularExpressionValidator()
        floatValidator.setRegularExpression(QRegularExpression("[A-Za-z]+"))

        self.seqText.setValidator(floatValidator)

        # length setting
        self.lenLabeLen1 = QLabel("1:")
        self.lenTextLen1 = QSpinBox()
        self.lenTextLen1.setMinimum(5)
        self.lenTextLen1.setMaximum(self.seqLen)
        self.skipLabeLen1 = QLabel("Skip:")
        self.skipTextLen1 = QSpinBox()
        self.skipTextLen1.setMinimum(0)
        self.skipTextLen1.setMaximum(self.seqLen)
        self.tailLabeLen1 = QLabel("Tail-off:")
        self.tailTextLen1 = QSpinBox()
        self.tailTextLen1.setMinimum(0)
        self.tailTextLen1.setMaximum(self.seqLen)
        self.offsetLabeLen1 = QLabel("Offset:")
        self.offsetTextLen1 = QSpinBox()
        self.offsetTextLen1.setMinimum(1)
        self.offsetTextLen1.setMaximum(10)

        self.lenLabeLen2 = QLabel("  2:")
        self.lenTextLen2 = QSpinBox()
        self.lenTextLen2.setMinimum(5)
        self.lenTextLen2.setMaximum(self.seqLen)
        self.skipLabeLen2 = QLabel("Skip:")
        self.skipTextLen2 = QSpinBox()
        self.skipTextLen2.setMinimum(0)
        self.skipTextLen2.setMaximum(self.seqLen)
        self.tailLabeLen2 = QLabel("Tail-off:")
        self.tailTextLen2 = QSpinBox()
        self.tailTextLen2.setMinimum(0)
        self.tailTextLen2.setMaximum(self.seqLen)
        self.offsetLabeLen2 = QLabel("Offset:")
        self.offsetTextLen2 = QSpinBox()
        self.offsetTextLen2.setMinimum(1)
        self.offsetTextLen2.setMaximum(10)

        self.lenLabeLen3 = QLabel("3:")
        self.lenTextLen3 = QSpinBox()
        self.lenTextLen3.setMinimum(5)
        self.lenTextLen3.setMaximum(35)
        self.skipLabeLen3 = QLabel("Skip:")
        self.skipTextLen3 = QSpinBox()
        self.skipTextLen3.setMinimum(0)
        self.skipTextLen3.setMaximum(100)
        self.tailLabeLen3 = QLabel("Tail-off:")
        self.tailTextLen3 = QSpinBox()
        self.tailTextLen3.setMinimum(0)
        self.tailTextLen3.setMaximum(100)
        self.offsetLabeLen3 = QLabel("Offset:")
        self.offsetTextLen3 = QSpinBox()
        self.offsetTextLen3.setMinimum(1)
        self.offsetTextLen3.setMaximum(10)

        self.lenLabeLen4 = QLabel("  4:")
        self.lenTextLen4 = QSpinBox()
        self.lenTextLen4.setMinimum(5)
        self.lenTextLen4.setMaximum(35)
        self.skipLabeLen4 = QLabel("Skip:")
        self.skipTextLen4 = QSpinBox()
        self.skipTextLen4.setMinimum(0)
        self.skipTextLen4.setMaximum(100)
        self.tailLabeLen4 = QLabel("Tail-off:")
        self.tailTextLen4 = QSpinBox()
        self.tailTextLen4.setMinimum(0)
        self.tailTextLen4.setMaximum(100)
        self.offsetLabeLen4 = QLabel("Offset:")
        self.offsetTextLen4 = QSpinBox()
        self.offsetTextLen4.setMinimum(1)
        self.offsetTextLen4.setMaximum(10)

        self.lenLabeLen5 = QLabel("5:")
        self.lenTextLen5 = QSpinBox()
        self.lenTextLen5.setMinimum(5)
        self.lenTextLen5.setMaximum(35)
        self.skipLabeLen5 = QLabel("Skip:")
        self.skipTextLen5 = QSpinBox()
        self.skipTextLen5.setMinimum(0)
        self.skipTextLen5.setMaximum(100)
        self.tailLabeLen5 = QLabel("Tail-off:")
        self.tailTextLen5 = QSpinBox()
        self.tailTextLen5.setMinimum(0)
        self.tailTextLen5.setMaximum(100)
        self.offsetLabeLen5 = QLabel("Offset:")
        self.offsetTextLen5 = QSpinBox()
        self.offsetTextLen5.setMinimum(1)
        self.offsetTextLen5.setMaximum(10)

        self.lenLabeLen6 = QLabel("  6:")
        self.lenTextLen6 = QSpinBox()
        self.lenTextLen6.setMinimum(5)
        self.lenTextLen6.setMaximum(35)
        self.skipLabeLen6 = QLabel("Skip:")
        self.skipTextLen6 = QSpinBox()
        self.skipTextLen6.setMinimum(0)
        self.skipTextLen6.setMaximum(100)
        self.tailLabeLen6 = QLabel("Tail-off:")
        self.tailTextLen6 = QSpinBox()
        self.tailTextLen6.setMinimum(0)
        self.tailTextLen6.setMaximum(100)
        self.offsetLabeLen6 = QLabel("Offset:")
        self.offsetTextLen6 = QSpinBox()
        self.offsetTextLen6.setMinimum(1)
        self.offsetTextLen6.setMaximum(10)

        self.peptypeInnerRandomizeCb = QCheckBox('Randomize pep for each peptype')
        self.peptypeInnerRandomizeCb.setChecked(True)
        self.btSubmit = QPushButton("SUBMIT")
        self.btSubmit.clicked.connect(self.submitBtClicked)
        self.btClear = QPushButton('CLEAR')
        self.btClear.clicked.connect(self.clearBtClicked)

        zeroRowLayout.addWidget(self.userLabel, 0,0,1,1,Qt.AlignmentFlag.AlignRight)
        zeroRowLayout.addWidget(self.userCombo, 0,1,1,1)
        zeroRowLayout.addWidget(self.projNameLabel, 0,2,1,1,Qt.AlignmentFlag.AlignRight)
        zeroRowLayout.addWidget(self.projNameText, 0,3,1,1)
        zeroRowLayout.addWidget(self.synNumLabel, 0,4,1,1,Qt.AlignmentFlag.AlignRight)
        zeroRowLayout.addWidget(self.synNumText, 0,5,1,1)
        zeroRowLayout.addWidget(self.synTypeLabel, 0, 6, 1, 1, Qt.AlignmentFlag.AlignRight)
        zeroRowLayout.addWidget(self.synTypeText, 0, 7, 1, 1)


        zeroRowLayout.addWidget(self.targetLabel, 1, 0, 1, 1,Qt.AlignmentFlag.AlignRight)
        zeroRowLayout.addWidget(self.targetText, 1, 1, 1, 1)
        zeroRowLayout.addWidget(self.pdbLabel, 1, 2, 1, 1,Qt.AlignmentFlag.AlignRight)
        zeroRowLayout.addWidget(self.pdbText, 1, 3, 1, 1)
        zeroRowLayout.addWidget(self.nSeqLabel, 1, 4, 1, 1,Qt.AlignmentFlag.AlignRight)
        zeroRowLayout.addWidget(self.nSeqText, 1, 5, 1, 1)
        zeroRowLayout.addWidget(self.seqNrLabel, 1, 6, 1, 1, Qt.AlignmentFlag.AlignRight)
        zeroRowLayout.addWidget(self.seqNrText, 1, 7, 1, 1)

        thirdRowLayout.addWidget(self.seqLabel)
        thirdRowLayout.addWidget(self.seqText,10)
        thirdRowLayout.addWidget(self.seqLenText)

        # length layout
        lenLayout.addWidget(self.lenLabeLen1,0,0,Qt.AlignmentFlag.AlignRight)
        lenLayout.addWidget(self.lenTextLen1,0,1)
        lenLayout.addWidget(self.skipLabeLen1,0,2)
        lenLayout.addWidget(self.skipTextLen1,0,3)
        lenLayout.addWidget(self.tailLabeLen1,0,4)
        lenLayout.addWidget(self.tailTextLen1,0,5)
        lenLayout.addWidget(self.offsetLabeLen1,0,6)
        lenLayout.addWidget(self.offsetTextLen1,0,7)
        #
        self.separatorLine = QFrame()
        self.separatorLine.setFrameShape(QFrame.Shape.VLine)
        lenLayout.addWidget(self.separatorLine,0,8)

        lenLayout.addWidget(self.lenLabeLen2,0,8,Qt.AlignmentFlag.AlignRight)
        lenLayout.addWidget(self.lenTextLen2,0,9)
        lenLayout.addWidget(self.skipLabeLen2,0,10)
        lenLayout.addWidget(self.skipTextLen2,0,11)
        lenLayout.addWidget(self.tailLabeLen2,0,12)
        lenLayout.addWidget(self.tailTextLen2,0,13)
        lenLayout.addWidget(self.offsetLabeLen2,0,14)
        lenLayout.addWidget(self.offsetTextLen2,0,15)

        lenLayout.addWidget(self.lenLabeLen3, 1, 0, Qt.AlignmentFlag.AlignRight)
        lenLayout.addWidget(self.lenTextLen3, 1, 1)
        lenLayout.addWidget(self.skipLabeLen3, 1, 2)
        lenLayout.addWidget(self.skipTextLen3, 1, 3)
        lenLayout.addWidget(self.tailLabeLen3,1,4)
        lenLayout.addWidget(self.tailTextLen3,1,5)
        lenLayout.addWidget(self.offsetLabeLen3, 1, 6)
        lenLayout.addWidget(self.offsetTextLen3, 1, 7)

        self.separatorLine = QFrame()
        self.separatorLine.setFrameShape(QFrame.Shape.VLine)
        lenLayout.addWidget(self.separatorLine, 1, 8)

        lenLayout.addWidget(self.lenLabeLen4, 1, 8, Qt.AlignmentFlag.AlignRight)
        lenLayout.addWidget(self.lenTextLen4, 1, 9)
        lenLayout.addWidget(self.skipLabeLen4, 1, 10)
        lenLayout.addWidget(self.skipTextLen4, 1, 11)
        lenLayout.addWidget(self.tailLabeLen4,1,12)
        lenLayout.addWidget(self.tailTextLen4,1,13)
        lenLayout.addWidget(self.offsetLabeLen4, 1, 14)
        lenLayout.addWidget(self.offsetTextLen4, 1, 15)

        lenLayout.addWidget(self.lenLabeLen5, 2, 0, Qt.AlignmentFlag.AlignRight)
        lenLayout.addWidget(self.lenTextLen5, 2, 1)
        lenLayout.addWidget(self.skipLabeLen5, 2, 2)
        lenLayout.addWidget(self.skipTextLen5, 2, 3)
        lenLayout.addWidget(self.tailLabeLen5,2,4)
        lenLayout.addWidget(self.tailTextLen5,2,5)
        lenLayout.addWidget(self.offsetLabeLen5, 2, 6)
        lenLayout.addWidget(self.offsetTextLen5, 2, 7)

        self.separatorLine = QFrame()
        self.separatorLine.setFrameShape(QFrame.Shape.VLine)
        lenLayout.addWidget(self.separatorLine,2,8)

        lenLayout.addWidget(self.lenLabeLen6, 2, 8, Qt.AlignmentFlag.AlignRight)
        lenLayout.addWidget(self.lenTextLen6, 2, 9)
        lenLayout.addWidget(self.skipLabeLen6, 2, 10)
        lenLayout.addWidget(self.skipTextLen6, 2, 11)
        lenLayout.addWidget(self.tailLabeLen6,2,12)
        lenLayout.addWidget(self.tailTextLen6,2,13)
        lenLayout.addWidget(self.offsetLabeLen6, 2, 14)
        lenLayout.addWidget(self.offsetTextLen6, 2, 15)

        self.browser =QWebEngineView()
        self.html = """
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            """
        self.browser.setHtml(self.html)
        self.browser.setFixedHeight(330)
        pdbLayout.addWidget(self.browser)

        btLayout.addStretch()
        btLayout.addWidget(self.btSubmit)
        btLayout.addWidget(self.peptypeInnerRandomizeCb)
        btLayout.addStretch()
        btLayout.addWidget(self.btClear)
        btLayout.addStretch()

        self.lenTextList = [self.lenTextLen1,self.lenTextLen2,self.lenTextLen3,self.lenTextLen4,
                       self.lenTextLen5,self.lenTextLen6]
        self.skipTextList = [self.skipTextLen1,self.skipTextLen2,self.skipTextLen3,self.skipTextLen4,
                       self.skipTextLen5,self.skipTextLen6]
        self.tailTextList = [self.tailTextLen1, self.tailTextLen2, self.tailTextLen3, self.tailTextLen4,
                             self.tailTextLen5, self.tailTextLen6]
        self.offsetTextList = [self.offsetTextLen1, self.offsetTextLen2, self.offsetTextLen3, self.offsetTextLen4,
                        self.offsetTextLen5, self.offsetTextLen6]

        for i in range(1,6):
            self.lenTextList[i].setEnabled(False)
            self.skipTextList[i].setEnabled(False)
            self.tailTextList[i].setEnabled(False)
            self.offsetTextList[i].setEnabled(False)

        layoutStack1.addWidget(zeroRow_Group)
        layoutStack1.addWidget(pdb_Group)
        layoutStack1.addLayout(thirdRowLayout)
        layoutStack1.addWidget(lenGroup)
        layoutStack1.addLayout(btLayout)

        self.stack1.setLayout(layoutStack1)

    def nLengthChanged(self,sig):
        for i in range(6):
            self.lenTextList[i].setEnabled(False)
            self.skipTextList[i].setEnabled(False)
            self.tailTextList[i].setEnabled(False)
            self.offsetTextList[i].setEnabled(False)
        for i in range(sig):
            self.lenTextList[i].setEnabled(True)
            self.skipTextList[i].setEnabled(True)
            self.tailTextList[i].setEnabled(True)
            self.offsetTextList[i].setEnabled(True)

    def seqPasted(self):
        self.seqLen = len(self.seqText.text())
        self.seqLenText.setText(str(self.seqLen))
        for i in range(6):
            self.skipTextList[i].setMaximum(self.seqLen)
            self.lenTextList[i].setMaximum(self.seqLen)
            self.skipTextList[i].valueChanged.connect(self.skipChanged)

    def skipChanged(self,v):
        mybox = self.sender()
        for i in range(6):
            if mybox is self.skipTextList[i]:
                self.tailTextList[i].setMaximum(self.seqLen - v)
                self.offsetTextList[i].setMaximum(self.seqLen - self.skipTextList[i].value())

    
    def pdbTextChanged(self):
        text = self.pdbText.text()
        self.html = """
                <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
                    """
        request = requests.get(f"https://files.rcsb.org/view/{text}.pdb")
        responseCode = request.status_code
        if responseCode == 404:
            self.html += "<h2 style='color:purple'> No structure available in PDB database </h2>"
        self.html += """ 
        <div style="height: 800px; width: 800px; position: relative;" class='viewer_3Dmoljs' 
        data-backgroundcolor='0xffffff' data-select1='chain:A' data-style1='cartoon:color=spectrum' 
        data-surface1='opacity:.5;color:white' data-select2='chain:B' data-style2='stick'
        """

        self.html += f"data-pdb='{text}'"
        self.html += "> </div> "
        self.browser.setHtml(self.html)

    
    def clearBtClicked(self):
        self.projNameText.clear() 
        self.targetText.clear()  
        self.pdbText.clear()   
        self.seqText.clear()
        self.nSeqText.setValue(1)
        self.synNumText.clear()
        self.browser.setHtml(" ")

    def generatePeptide(self, pepSeq=" ", pepLength=None, offset=None, tailoff = None, skip=None):
        sequence = pepSeq[int(skip):(len(pepSeq)-int(tailoff))]
        sequence = sequence.upper()
        pepSequenceList = []
        i = 0
        while i in range(len(sequence) - int(pepLength) + 1):
            pepSequence = sequence[i:(i + int(pepLength))]
            pepSequenceList.append(pepSequence)
            i += int(offset)
        pepSequenceTable = pd.DataFrame(pepSequenceList,columns=['SEQUENCE'])
        # np.random.seed(108)
        # if randomize == True:
        #     pepSequenceTable = pepSequenceTable.sample(frac=1, random_state=108)
        return pepSequenceTable

    def submitBtClicked(self):
        self.user = self.userCombo.currentText()
        self.projectName = self.projNameText.text()
        self.synthesisNumber = self.synNumText.text()
        self.synType = self.synTypeText.currentText()
        self.target = self.targetText.text()
        self.nLengths = int(self.nSeqText.text())
        self.sequence = self.seqText.text()

        paraList = [self.user,self.projectName,self.synthesisNumber,self.synType,self.target,self.sequence]

        if '' in paraList:
            QMessageBox.information(self, 'Info', 'Please set all parameters', QMessageBox.StandardButton.Ok)
        else:
            self.lenListFinal = []
            self.skipListFinal = []
            self.tailListFinal = []
            self.offsetListFinal = []
            self.typenameFinal = []
            for i in range(self.nLengths):
                name = self.synType + str(self.lenTextList[i].value())
                self.typenameFinal.append(name)
                self.lenListFinal.append(self.lenTextList[i].value())
                self.skipListFinal.append(self.skipTextList[i].value())
                self.tailListFinal.append(self.tailTextList[i].value())
                self.offsetListFinal.append(self.offsetTextList[i].value())
            self.dfLen = pd.concat([pd.Series(self.lenListFinal), pd.Series(self.skipListFinal),pd.Series(self.tailListFinal),pd.Series(self.offsetListFinal)], axis=1)
            self.dfLen.index = self.typenameFinal
            self.dfLen.columns = ['Length','Skip','Tailoff','Offset']
            self.sequenceTableStackTemp = []
            self.overviewDict= dict()
            for i in range(len(self.typenameFinal)):
                sequenceTable = self.generatePeptide(pepSeq=self.sequence, pepLength=self.lenListFinal[i],skip=self.skipListFinal[i],
                                                    offset=self.offsetListFinal[i],tailoff=self.tailListFinal[i])
                sequenceTable['infoName'] = [(str(item)+f"-{self.typenameFinal[i]}-offset{self.offsetListFinal[i]}" ) for item in sequenceTable.index]
                sequenceTable['peptide_label'] = self.typenameFinal[i]
                sequenceTable['Syn'] = self.synthesisNumber
                sequenceTable['project_naam'] = self.projectName
                sequenceTable['contact_persoon'] = self.user
                sequenceTable['groep_type'] = self.typenameFinal[i]
                sequenceTable['target_protein'] = self.target
                sequenceTable['n_capping'] = 'Dummy material'
                sequenceTable['acrylic_acid'] = '88%'
                sequenceTable['graft_batch'] = 40
                sequenceTable['clips_type'] = self.typenameFinal[i]
                sequenceTable['Brief'] = self.projectName
                sequenceTable['Label'] = self.typenameFinal[i]
                sequenceTable['ProjectID'] = self.projectName
                sequenceTable['ContactPersoon'] = self.user
                sequenceTable['AcrylicAcid'] = '88%'
                sequenceTable['CLIPS'] = self.typenameFinal[i]
                sequenceTable['n.capping'] = 'Dummy material2'
                if self.peptypeInnerRandomizeCb.checkState() == Qt.CheckState.Checked:
                    np.random.seed(108)
                    sequenceTable = sequenceTable.sample(frac=1, random_state=108)

                self.sequenceTableStackTemp.append(sequenceTable)
                overviewDictKey = self.typenameFinal[i]+'-offset'+str(self.offsetListFinal[i])
                self.overviewDict[overviewDictKey] = sequenceTable.shape[0]
            self.sequenceTableStack = pd.concat([item for item in self.sequenceTableStackTemp], axis=0, sort=False, ignore_index=True)
            self.NrSeq = self.seqNrText.value()

            self.pepToStoreVal = sum(self.overviewDict.values())
            self.cut1Nr.setMaximum(self.pepToStoreVal)
            self.cut1Nr.setValue(self.pepToStoreVal)
            self.cut2Nr.setMaximum(self.pepToStoreVal)
            self.cut2Nr.setValue(self.pepToStoreVal)


            self.sequenceTableStack['infoName'] = self.sequenceTableStack['infoName'] + f'-seq{self.NrSeq}'
            self.overviewDf = pd.DataFrame.from_dict(self.overviewDict,orient='index')
            self.overviewDf.columns=['Value']
            self.overviewDf.loc['Number_peptide'] = self.pepToStoreVal
            self.overviewDf.loc['Number_minicard'] = np.ceil(self.pepToStoreVal/455)
            self.overviewDf.loc['Number_TYRcontrol'] = 0
            self.overviewDf.loc['Number_PEPcontrol'] = 0
            self.seqTableStackKey = f'seq{self.NrSeq}-' + self.synType + f'-total{self.pepToStoreVal}'

            # else:
            #     QMessageBox.information(self, 'Info', 'Please input correct information', QMessageBox.StandardButton.Ok)

            self.overviewModel = TableModel(self.overviewDf)
            self.overviewTable.setModel(self.overviewModel)
            self.overviewTable.setColumnWidth(0, 400)
            self.list.setCurrentRow(1)

    def stack2Layout(self):
        self.layoutStack2 = QHBoxLayout()
        self.layoutStack2Sub = QVBoxLayout()
        self.layoutStack2Left = QVBoxLayout()
        layoutStack2Right = QVBoxLayout()
        controlLayout = QFormLayout()
        self.listsLayout = QHBoxLayout()
        leftListLayout = QVBoxLayout()
        rightListLayout = QVBoxLayout()

        self.overviewTable = QTableView()
        self.overviewTable.setColumnWidth(0,400)

        # self.sc = MplCanvas(self.layoutStack2Left, dpi=100)
        # toolbar = NavigationToolbar(self.sc, self)
        self.listLeft = QListWidget()
        self.listLeft.setDragEnabled(True)
        self.listLeft.setAcceptDrops(True)
        self.listLeft.setDefaultDropAction(Qt.DropAction.LinkAction)

        self.listRight = QListWidget()
        self.listRight.setAcceptDrops(True)
        self.listRight.setDragDropMode(QAbstractItemView.DragDropMode.DragDrop)
        self.listRight.setDefaultDropAction(Qt.DropAction.MoveAction)
        self.listRight.currentItemChanged.connect(self.listRightChanged)
        self.listRight.itemChanged.connect(self.listRightChanged)


        leftList_Group = QGroupBox(self)
        leftList_Group.setTitle("&Storehouse")
        leftList_Group.setLayout(leftListLayout)
        rightList_Group = QGroupBox(self)
        rightList_Group.setTitle("&Factory")
        rightList_Group.setLayout(rightListLayout)

        self.addToStoreHouseBt = QPushButton('Add peptide set to storehouse  ')
        self.addToStoreHouseBt.clicked.connect(self.addToStoreHouseBtClicked)
        self.addSeqBt = QPushButton('Add another peptype or sequence')
        self.addSeqBt.clicked.connect(self.addSeqBtClicked)


        cut1Lb = QLabel('cutPoint1')
        cut2Lb = QLabel('cutPoint2')
        self.cut1Nr = QSpinBox()
        self.cut1Nr.setMinimum(1)
        self.cut2Nr = QSpinBox()
        self.cut2Nr.setMinimum(1)
        self.innerRandomizeCb = QCheckBox('Mix peptype')
        self.innerRandomizeCb.setChecked(False)
        nrEnRandLayout = QHBoxLayout()
        nrEnRandLayout.addWidget(cut1Lb)
        nrEnRandLayout.addWidget(self.cut1Nr)
        nrEnRandLayout.addWidget(cut2Lb)
        nrEnRandLayout.addWidget(self.cut2Nr)

        self.generateCardBt = QPushButton('Generate minicards')
        self.generateCardBt.clicked.connect(self.generateCardBtClicked)
        self.tyrText = QSpinBox()
        self.tyrText.setValue(1)
        self.tyrText.setMinimum(1)
        self.tyrText.setMaximum(30)
        self.tyrBt = QPushButton('Add TYR control to storehouse')
        self.tyrBt.clicked.connect(self.addTYRBtClicked)
        self.pepCtlText = QSpinBox()
        self.pepCtlText.setValue(1)
        self.pepCtlText.setMinimum(1)
        self.pepCtlText.setMaximum(400)
        self.pepCtlBt = QPushButton('Add Pep control to storehouse')
        self.pepCtlBt.clicked.connect(self.pepCtlBtClicked)

        nrMiniCardsLb = QPushButton('Number of minicard in factory: ')
        self.nrMiniCardsText = QLabel()
        self.nrMiniCardsText.setText('0')
        self.nrMiniCardsText.setAlignment(Qt.AlignmentFlag.AlignCenter)

        nrPepLb = QPushButton('Number of peptides in factory: ')
        self.nrPepText = QLabel()
        self.nrPepText.setText('0')
        self.nrPepText.setAlignment(Qt.AlignmentFlag.AlignCenter)

        nrEmptyLb = QPushButton('    Wells available in last card: ')
        self.nrEmptyText = QLabel()
        self.nrEmptyText.setText('455')
        self.nrEmptyText.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.cycleLb = QLabel()
        self.cycleLb.setPixmap(QPixmap('cycle.png'))
        self.cycleLb.setScaledContents(True)
        self.cycleLb.setMaximumWidth(50)
        self.cycleLb.setMaximumHeight(50)

        self.clearLeftListBt = QPushButton('Clear')
        self.clearLeftListBt.clicked.connect(self.clearLeftListBtClicked)

        self.layoutStack2Left.addWidget(self.overviewTable,5)
        controlLayout.addRow(self.addToStoreHouseBt,nrEnRandLayout)
        controlLayout.addRow('',self.innerRandomizeCb)
        controlLayout.addRow(self.tyrBt,self.tyrText)
        controlLayout.addRow(self.pepCtlBt,self.pepCtlText)
        controlLayout.addRow(nrMiniCardsLb,self.nrMiniCardsText)
        controlLayout.addRow(nrPepLb,self.nrPepText)
        controlLayout.addRow(nrEmptyLb,self.nrEmptyText)

        # layoutStack2Right.addWidget(self.addToStoreHouseBt)
        layoutStack2Right.addLayout(controlLayout)
        layoutStack2Right.addWidget(self.addSeqBt)
        # layoutStack2Right.addWidget(self.generateCardBt)

        leftListLayout.addWidget(self.listLeft)
        leftListLayout.addWidget(self.clearLeftListBt)
        leftList_Group.setLayout(leftListLayout)
        rightListLayout.addWidget(self.listRight)
        rightListLayout.addWidget(self.generateCardBt)
        rightList_Group.setLayout(rightListLayout)

        self.listsLayout.addWidget(leftList_Group)
        self.listsLayout.addWidget(self.cycleLb)
        self.listsLayout.addWidget(rightList_Group)

        self.layoutStack2.addLayout(self.layoutStack2Left,8)
        self.layoutStack2.addStretch(1)
        self.layoutStack2.addLayout(layoutStack2Right,8)
        self.layoutStack2Sub.addLayout(self.layoutStack2,4)
        self.layoutStack2Sub.addLayout(self.listsLayout,6)
        self.stack2.setLayout(self.layoutStack2Sub)
        self.overviewTable.update()

        self.pepStoreDict = dict()

    def addToStoreHouseBtClicked(self):
        if hasattr(self, 'sequenceTableStack'):
            if self.innerRandomizeCb.checkState() == Qt.CheckState.Checked:
                np.random.seed(108)
                self.sequenceTableStack = self.sequenceTableStack.sample(frac=1, random_state=108)

            self.pepStoreDict[self.seqTableStackKey] = self.sequenceTableStack
            if self.cut1Nr.value() < self.pepToStoreVal and self.cut2Nr.value() == self.pepToStoreVal:
                cutNrBeforeToStore = self.cut1Nr.value()
                cutNrBeforeToStore2 = self.pepToStoreVal - self.cut1Nr.value()
                keyTempPart1 = self.seqTableStackKey + '-part1' + f'-{cutNrBeforeToStore}'
                keyTempPart2 = self.seqTableStackKey + '-part2' + f'-{cutNrBeforeToStore2}'
                self.pepStoreDict[keyTempPart1] = self.sequenceTableStack[:cutNrBeforeToStore]
                self.pepStoreDict[keyTempPart2] = self.sequenceTableStack[cutNrBeforeToStore:]
                self.listLeft.addItem(keyTempPart1)
                self.listLeft.addItem(keyTempPart2)
            elif self.cut1Nr.value() == self.pepToStoreVal and self.cut2Nr.value() < self.pepToStoreVal:
                QMessageBox.information(self, 'Info', 'Cutting into two parts? use cut1', QMessageBox.StandardButton.Ok)
                self.cut2Nr.setValue(self.pepToStoreVal)
            elif self.cut1Nr.value() < self.pepToStoreVal and self.cut2Nr.value() < self.pepToStoreVal and self.cut1Nr.value() < self.cut2Nr.value():
                cutNrBeforeToStore = self.cut1Nr.value()
                cutNrBeforeToStore2 = self.pepToStoreVal - self.cut1Nr.value()
                cutNrBeforeToStore3 = self.cut2Nr.value() - self.cut1Nr.value()
                cutNrBeforeToStore4 = self.pepToStoreVal - self.cut2Nr.value()
                keyTempPart1 = self.seqTableStackKey + '-part1' + f'-{cutNrBeforeToStore}'
                keyTempPart2 = self.seqTableStackKey + '-part2' + f'-{cutNrBeforeToStore3}'
                keyTempPart3 = self.seqTableStackKey + '-part3' + f'-{cutNrBeforeToStore4}'
                self.pepStoreDict[keyTempPart1] = self.sequenceTableStack[:cutNrBeforeToStore]
                self.pepStoreDict[keyTempPart2] = self.sequenceTableStack[cutNrBeforeToStore:self.cut2Nr.value()]
                self.pepStoreDict[keyTempPart3] = self.sequenceTableStack[self.cut2Nr.value():]
                self.listLeft.addItem(keyTempPart1)
                self.listLeft.addItem(keyTempPart2)
                self.listLeft.addItem(keyTempPart3)
            elif self.cut1Nr.value() < self.pepToStoreVal and self.cut2Nr.value() < self.pepToStoreVal and self.cut1Nr.value() >= self.cut2Nr.value():
                QMessageBox.information(self, 'Info', 'cut1 cannot be higher than cut2', QMessageBox.StandardButton.Ok)
            else:
                if self.listLeft.count() > 0 and (self.seqTableStackKey in [self.listLeft.item(i).text() for i in range(self.listLeft.count())]):
                    QMessageBox.information(self, 'Info', 'The selected sample is already in list', QMessageBox.StandardButton.Ok)
                else:
                    self.listLeft.addItem(self.seqTableStackKey)  
        else:
            QMessageBox.information(self, 'Info', 'Please generate peptides first', QMessageBox.StandardButton.Ok)

    def addTYRBtClicked(self):
        try:
            self.tyrNr = self.tyrText.value()
            tyrSeq = []
            for i in range(self.tyrNr):
                item = ['1'+'X'*i]
                tyrSeq.append(item)
            self.tyrDf = pd.DataFrame(tyrSeq)
            self.tyrDf.columns = ['SEQUENCE']
            self.tyrDf['infoName'] = [(str(item) + '-' + 'TYR'+'-'+'offset0'+'-'+f'seq{self.tyrNr}') for item in self.tyrDf.index]
            self.tyrDf['peptide_label'] = 'TYR'
            self.tyrDf['Syn'] = self.synthesisNumber
            self.tyrDf['project_naam'] = 'NewLinearControls'
            self.tyrDf['contact_persoon'] = self.user
            self.tyrDf['groep_type'] = 'Linear'
            self.tyrDf['target_protein'] = 'new Ctrls'
            self.tyrDf['n_capping'] = 'Dummy material3'
            self.tyrDf['acrylic_acid'] = '88%'
            self.tyrDf['graft_batch'] = 40
            self.tyrDf['clips_type'] = 'Linear'
            self.tyrDf['Brief'] = 'NewT3'
            self.tyrDf['Label'] = 'TYR'
            self.tyrDf['ProjectID'] = 'NewLinearControls'
            self.tyrDf['ContactPersoon'] = self.user
            self.tyrDf['AcrylicAcid'] = '88%'
            self.tyrDf['CLIPS'] = 'Linear'
            self.tyrDf['n.capping'] = 'Dummy material4'

            self.tyrKey = f'seq0-TYR-total{self.tyrNr}'
            self.pepStoreDict[self.tyrKey] = self.tyrDf
            if self.listLeft.count() > 0 and (self.tyrKey in [self.listLeft.item(i).text() for i in range(self.listLeft.count())]):
                QMessageBox.information(self, 'Info', 'The selected sample is already in list', QMessageBox.StandardButton.Ok)
            else:
                self.listLeft.addItem(self.tyrKey)
        except:
            QMessageBox.information(self, 'Info', 'Please fill in all parameters first', QMessageBox.StandardButton.Ok)

    def pepCtlBtClicked(self):
        try:
            self.pepCtlNr = self.pepCtlText.value()
            n = int(np.ceil(self.pepCtlNr/len(self.pepCtl)))
            self.pepCtlRep = self.pepCtl*n
            self.pepCtlSelected = self.pepCtlRep[:self.pepCtlNr]
            self.pepCtlDf = pd.DataFrame(self.pepCtlSelected)
            self.pepCtlDf.columns = ['SEQUENCE']
            self.pepCtlDf['infoName'] = [(str(item) + '-' + 'PEP' + '-' + 'offset0' + '-' + f'seq{self.pepCtlNr}') for item in
                                    self.pepCtlDf.index]
            self.pepCtlDf['peptide_label'] = 'PEP'
            self.pepCtlDf['Syn'] = self.synthesisNumber
            self.pepCtlDf['project_naam'] = 'NewLinearControls'
            self.pepCtlDf['contact_persoon'] = self.user
            self.pepCtlDf['groep_type'] = 'Linear'
            self.pepCtlDf['target_protein'] = 'new Ctrls'
            self.pepCtlDf['n_capping'] = 'Dummy material3'
            self.pepCtlDf['acrylic_acid'] = '88%'
            self.pepCtlDf['graft_batch'] = 40
            self.pepCtlDf['clips_type'] = 'Linear'
            self.pepCtlDf['Brief'] = 'NewT3'
            self.pepCtlDf['Label'] = 'TYR'
            self.pepCtlDf['ProjectID'] = 'NewLinearControls'
            self.pepCtlDf['ContactPersoon'] = self.user
            self.pepCtlDf['AcrylicAcid'] = '88%'
            self.pepCtlDf['CLIPS'] = 'Linear'
            self.pepCtlDf['n.capping'] = 'Dummy material4'

            self.pepKey = f'seq0-PEP-total{self.pepCtlNr}'
            self.pepStoreDict[self.pepKey] = self.pepCtlDf
            if self.listLeft.count() > 0 and (self.pepKey in [self.listLeft.item(i).text() for i in range(self.listLeft.count())]):
                QMessageBox.information(self, 'Info', 'The selected control is already in list', QMessageBox.StandardButton.Ok)
            else:
                self.listLeft.addItem(self.pepKey)
        except:
            QMessageBox.information(self, 'Info', 'Please fill in all parameters first', QMessageBox.StandardButton.Ok)

    def addSeqBtClicked(self):
        self.list.setCurrentRow(0)
        self.clearBtClicked()

    def clearLeftListBtClicked(self):
        self.listLeft.clear()

    def listRightChanged(self):
        itemList = [self.listRight.item(i).text() for i in range(self.listRight.count())]
        self.pepInFacNr = sum([self.pepStoreDict[item].shape[0] for item in itemList])
        self.nrPepText.setText(str(self.pepInFacNr))
        self.miniCardNr = np.ceil(self.pepInFacNr / 455)
        self.nrMiniCardsText.setText(str(self.miniCardNr))
        availableWells = self.miniCardNr*455 - self.pepInFacNr
        self.nrEmptyText.setText(str(availableWells))

    def generateCardBtClicked(self):
        if self.listRight.count() > 0:
            itemList = [self.listRight.item(i).text() for i in range(self.listRight.count())]
            dataList = []
            self.miniCardNr = int(self.miniCardNr)
            self.pepInFacNr = sum([self.pepStoreDict[item].shape[0] for item in itemList])
            for item in itemList:
                dataList.append(self.pepStoreDict[item])
            self.finalDf = pd.concat(dataList)
            cardID = sum([[i]*455 for i in range(1,self.miniCardNr+1)],[])
            cardID = cardID[:self.pepInFacNr]
            cardsynpos = list(range(1,456))*self.miniCardNr
            cardsynpos = cardsynpos[:self.pepInFacNr]
            m = int(np.ceil(self.pepInFacNr/17))
            cardrow = list(range(1,18))*(m+1)
            del cardrow[::459]
            del cardrow[15::458]
            del cardrow[440::457]
            del cardrow[455::456]
            cardrow = cardrow[:self.pepInFacNr]
            cardcol = sum([[i]*17 for i in range(1,m+2)],[])
            del cardcol[::459]
            del cardcol[15::458]
            del cardcol[440::457]
            del cardcol[455::456]
            cardcol = cardcol[:self.pepInFacNr]
            self.finalDf['cardID'] = cardID
            self.finalDf['card.synpos'] = cardsynpos
            self.finalDf['card.row'] = cardrow
            self.finalDf['card.col'] = cardcol
            self.finalDf = self.finalDf[['Syn','cardID','card.synpos','card.row','card.col','peptide_label','project_naam',
                                        'contact_persoon','SEQUENCE','groep_type','target_protein','n_capping',
                                        'acrylic_acid','graft_batch','clips_type','Brief','Label','ProjectID',
                                        'ContactPersoon','AcrylicAcid','CLIPS','n.capping','infoName']]
            # if self.finalDf.shape[0] > 1:
            #     saveFileName = QFileDialog.getSaveFileName(self, 'Save file','','Excel files (*.xlsx);; CSV (*.csv)')
            #     if saveFileName and saveFileName[1] == 'Excel files (*.xlsx)':
            #         self.finalDf.to_excel(saveFileName[0])
            #     elif saveFileName and saveFileName[1] == 'CSV (*.csv)':
            #         self.finalDf.to_csv(saveFileName[0])
            #     else:
            #         pass
            # else:
            #     QMessageBox.warning(self, 'Info','No data file can be downloaded', QMessageBox.StandardButton.Ok)

            figRow = int(np.ceil(self.miniCardNr / 2))
            if figRow <2:
                figRow = 2
            z1 = [0, 1] * figRow
            z2 = sum([[i, i] for i in range(figRow)], [])
            axList = [*zip(z2, z1)]

            try:
                for i in range(4):
                    self.scStack3.ax[axList[i][0], axList[i][1]].clear()
            except:
                pass

            if hasattr(self,'cardDict'):
                lenCardDict = len(self.cardDict)
            else:
                self.cardDict = dict()
                lenCardDict = 0

            for i in range(self.miniCardNr):
                ax = self.scStack3.ax[axList[i][0],axList[i][1]]
                scatterplot(x=self.finalDf['card.col'][i*455:455*(i+1)]-27*i, y=self.finalDf['card.row'][i*455:455*(i+1)], ax=ax,
                            hue=self.finalDf['peptide_label'][i*455:455*(i+1)], s=60, palette='tab20')
                ax.legend(ncol=1,loc='right',bbox_to_anchor=(1.2,0.5),prop={'size': 8})
                ax.set_xlim(0, 28)
                ax.set_ylim(0, 18)
                ax.invert_yaxis()
                ax.invert_xaxis()
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_title(f'Card-{lenCardDict+i+1}')

                self.finalDf[i * 455:455 * (i + 1)] = self.finalDf[i * 455:455 * (i + 1)].assign(cardID = lenCardDict+i+1)
                key = f'Card-{lenCardDict+i+1}'
                value = self.finalDf.iloc[i*455:455*(i+1),:]
                self.cardDict[key] = value

            plt.tight_layout()
            self.scStack3.draw()
            self.list.setCurrentRow(2)
        else:
            QMessageBox.warning(self, 'Info','No data file to generate cards', QMessageBox.StandardButton.Ok)

    def stack3Layout(self):
        self.layoutStack3 = QVBoxLayout()
        layoutArrNr = QFormLayout()
        layoutStack3Bottom = QHBoxLayout()

        self.scStack3 = MplCanvas(self.layoutStack2Left, dpi=100, nrows=2, ncols=2)
        toolbar = NavigationToolbar(self.scStack3, self)
        self.addToArrayBt = QPushButton('Submit for array generation')
        self.addToArrayBt.clicked.connect(self.addToArrayBtClicked)

        self.arrayNrLb = QLabel('Number of array: ')
        self.arrayNr = QSpinBox()
        self.arrayNr.setMinimum(1)
        self.arrayNr.setMaximum(3)
        layoutArrNr.addRow(self.arrayNrLb,self.arrayNr)

        layoutStack3Bottom.addStretch()
        layoutStack3Bottom.addLayout(layoutArrNr)
        layoutStack3Bottom.addWidget(self.addToArrayBt)
        layoutStack3Bottom.addStretch()
        self.layoutStack3.addWidget(toolbar)
        self.layoutStack3.addWidget(self.scStack3)
        self.layoutStack3.addLayout(layoutStack3Bottom)
        self.layoutStack3.addSpacing(10)
        self.stack3.setLayout(self.layoutStack3)

    def addToArrayBtClicked(self):
        if hasattr(self, 'cardDict'):
            nrOfArr = self.arrayNr.value()
            if nrOfArr == 2:
                self.array2List.setEnabled(True)
                self.generateCardBtClicked()
                withNewKeys = list(self.cardDict.keys())
                for i in range(self.miniCardNr):
                    self.infoText = self.infoText + f'INFOTEXT {withNewKeys[-(i+1)]} = {withNewKeys[-(i+1+self.miniCardNr)]} \n'
            elif nrOfArr == 3:
                self.array2List.setEnabled(True)
                self.array3List.setEnabled(True)
                self.generateCardBtClicked()
                self.generateCardBtClicked()
                withNewKeys = list(self.cardDict.keys())
                for i in range(self.miniCardNr):
                    self.infoText = self.infoText + f'INFOTEXT {withNewKeys[-(i+1)]} = {withNewKeys[-(i+1+self.miniCardNr)]} = {withNewKeys[-(i+1+self.miniCardNr*2)]}\n'
            self.cardList.clear()
            self.cardList.addItems(self.cardDict.keys())
            self.dupText.setText(self.infoText)
        else:
            QMessageBox.warning(self, 'Info','No cards available for array generation', QMessageBox.StandardButton.Ok)

    def stack4Layout(self):
        layoutStack4 = QHBoxLayout()
        stack4LeftListLayout = QVBoxLayout()
        stack4RightListLayout = QVBoxLayout()
        stack4BtLayout = QVBoxLayout()

        self.dupText = QTextEdit()

        self.cardList = QListWidget()
        self.cardList.setAcceptDrops(True)
        self.cardList.setDragDropMode(QAbstractItemView.DragDropMode.DragDrop)
        self.cardList.setDefaultDropAction(Qt.DropAction.MoveAction)
        cardList_Group = QGroupBox(self)
        cardList_Group.setTitle("&MiniCards")
        stack4LeftListLayout.addWidget(self.cardList)
        stack4LeftListLayout.addWidget(self.dupText)
        cardList_Group.setLayout(stack4LeftListLayout)

        array1Lb = QLabel('ARRAY-ONE')
        array2Lb = QLabel('ARRAY-TWO')
        array3Lb = QLabel('ARRAY-THREE')

        self.array1List = QListWidget()
        self.array1List.setAcceptDrops(True)
        self.array1List.setDragDropMode(QAbstractItemView.DragDropMode.DragDrop)
        self.array1List.setDefaultDropAction(Qt.DropAction.MoveAction)
        self.array2List = QListWidget()
        self.array2List.setAcceptDrops(True)
        self.array2List.setDragDropMode(QAbstractItemView.DragDropMode.DragDrop)
        self.array2List.setDefaultDropAction(Qt.DropAction.MoveAction)
        self.array2List.setEnabled(False)
        self.array3List = QListWidget()
        self.array3List.setAcceptDrops(True)
        self.array3List.setDragDropMode(QAbstractItemView.DragDropMode.DragDrop)
        self.array3List.setDefaultDropAction(Qt.DropAction.MoveAction)
        self.array3List.setEnabled(False)
        stack4RightListLayout.addWidget(array1Lb)
        stack4RightListLayout.addWidget(self.array1List)
        stack4RightListLayout.addWidget(array2Lb)
        stack4RightListLayout.addWidget(self.array2List)
        stack4RightListLayout.addWidget(array3Lb)
        stack4RightListLayout.addWidget(self.array3List)
        right_Group = QGroupBox(self)
        right_Group.setTitle("&Arrays")
        right_Group.setLayout(stack4RightListLayout)

        self.downldFinalBt = QPushButton('Download Final file')
        self.downldFinalBt.clicked.connect(self.downldFinalBtClicked)
        stack4BtLayout.addStretch(1)
        stack4BtLayout.addWidget(self.downldFinalBt)
        stack4BtLayout.addStretch(9)

        layoutStack4.addWidget(cardList_Group)
        layoutStack4.addWidget(right_Group)
        layoutStack4.addLayout(stack4BtLayout)
        self.stack4.setLayout(layoutStack4)

    def downldFinalBtClicked(self):
        self.array1ListKey = [self.array1List.item(i).text() for i in range(self.array1List.count())]
        if len(self.array1ListKey) > 0:
            for item in self.array1ListKey:
                 self.cardDict[item] = self.cardDict[item].assign(infoName = self.cardDict[item]['infoName']+ '-arr1')
        
        self.array2ListKey = [self.array2List.item(i).text() for i in range(self.array2List.count())]
        if len(self.array2ListKey) > 0:
            for item in self.array2ListKey:
                 self.cardDict[item] = self.cardDict[item].assign(infoName = self.cardDict[item]['infoName']+ '-arr2')

        self.array3ListKey = [self.array3List.item(i).text() for i in range(self.array3List.count())]
        if len(self.array3ListKey) > 0:
            for item in self.array3ListKey:
                 self.cardDict[item] = self.cardDict[item].assign(infoName = self.cardDict[item]['infoName']+ '-arr3')

        totalList = self.array1ListKey + self.array2ListKey + self.array3ListKey
        finalDfTemp = []
        for item in totalList:
            finalDfTemp.append(self.cardDict[item])
        
        if len(finalDfTemp) == 0:
            QMessageBox.warning(self, 'Info','No data file can be downloaded', QMessageBox.StandardButton.Ok)
        else:
            self.finalDfFile = pd.concat(finalDfTemp)

        if hasattr(self, 'finalDfFile') and self.finalDfFile.shape[0] > 1:
            saveFileName = QFileDialog.getSaveFileName(self, 'Save file','','Excel files (*.xlsx);; CSV (*.csv)')
            if saveFileName and saveFileName[1] == 'Excel files (*.xlsx)':
                self.finalDfFile.to_excel(saveFileName[0])
            elif saveFileName and saveFileName[1] == 'CSV (*.csv)':
                self.finalDfFile.to_csv(saveFileName[0])
            else:
                pass
        else:
            QMessageBox.warning(self, 'Info','No data file can be downloaded', QMessageBox.StandardButton.Ok)

if __name__ == '__main__':
    # def exceptionHook(exc_type, exc_value, exc_tb):
    #     tb = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
    #     currentTime = datetime.datetime.today().time().strftime('%H-%M-%S')
    #     filename = f"../{currentTime}-pepsynther-logFile.log"
    #     with open(filename,'w') as f:
    #         f.write(f'Biosynth pepsynther traceback info: \n {tb}')
    #     QMessageBox.warning(window,'warning', 'Something went wrong, please check log file for details', QMessageBox.StandardButton.Ok)
    #
    # sys.excepthook = exceptionHook
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    apply_stylesheet(app, theme='dark_amber.xml')
    window = PepSyn()
    window.setWindowIcon(QIcon('images/icoLK.ico'))
    window.show()
    sys.exit(app.exec())