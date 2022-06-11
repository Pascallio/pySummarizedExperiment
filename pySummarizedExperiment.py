import pandas as pd
import copy
import numpy as np

class pySummarizedExperiment():
    def __init__(self, assays: dict = dict(), columnData: pd.DataFrame = pd.DataFrame(), 
                 rowData: pd.DataFrame = pd.DataFrame(), metadata: dict = dict(),
                 longDf: pd.DataFrame = pd.DataFrame(), rowIndex: str = None, colIndex: str = None):

        if len(longDf.index) > 0:
            rowData, columnData, assays, metadata = self.__longToExperiment(longDf, rowIndex, colIndex)
            
        self.__assayData = {str(key): assays[key] for key in assays.keys()}
        self.__colDat = columnData
        self.__rowDat = rowData
        self.__metadat = metadata
        
    def __longToExperiment(self, df, rowIndex, colIndex):
        metadata = np.unique(df[df.columns[np.where(df.nunique() == 1)]])
        df = df[df.columns[np.where(df.nunique() > 1)]]
        df[rowIndex] = df[rowIndex].astype(str)
        df[colIndex] = df[colIndex].astype(str)
        rowCols = [rowIndex]
        for col in df.columns:
            if col != rowIndex:
                comb = df[[rowIndex, col]].drop_duplicates()
                one = len(set(comb[rowIndex]))                
                if one == len(comb.index):
                    rowCols.append(col)
        rowData = df[rowCols].drop_duplicates().set_index(rowIndex)
        
                    
        colCols = [colIndex]
        for col in df.columns:
            if col != colIndex:
                comb = df[[colIndex, col]].drop_duplicates()
                one = len(set(comb[colIndex]))
                if one == len(comb.index):
                    colCols.append(col)
        colData = df[colCols].drop_duplicates().set_index(colIndex)

        assays = {name: df[[name, rowIndex, colIndex]].pivot(index = rowIndex, columns = colIndex, values = name) 
                  for name in df.columns if name not in rowCols + colCols}
        return rowData, colData, assays, metadata 
            
    def validateWith(self, exp, axis):
        if axis == 0:
            cols = all([ind in exp.colnames() for ind in self.colnames()])
            cols2 = all([ind in exp.colData().columns for ind in self.colData().columns])

            rows = all([ind not in exp.index() for ind in self.index()])
            return cols and cols2 and rows
        
        cols = all([ind not in exp.colnames() for ind in self.colnames()])
        rows = all([ind in exp.index() for ind in self.index()])
        rows2 = all([ind in exp.rowData().columns for ind in self.rowData().columns])

        return cols and rows and rows2
    
    def bindRowData(self, exp):
        self.__rowDat = pd.concat([self.rowData(), exp.rowData()])
        for key in self.__assayData.keys():
            self.__assayData[key] = pd.concat([self.assay(key), exp.assay(key)], axis = 0)        

    def bindColData(self, exp):
        self.__colDat =  pd.concat([self.colData(), exp.colData()])
        for key in self.__assayData.keys():
            self.__assayData[key] = pd.concat([self.assay(key), exp.assay(key)], axis = 1)    
        
    def bind(self, exp, axis = 0):
        if axis not in [0, 1]:
            raise Exception("axis must be in range [0, 1]")
        x = copy.deepcopy(self)
        if x.validateWith(exp, axis):
            if axis == 0:
                x.bindRowData(exp)
            else: 
                x.bindColData(exp)
        return x

        
    def __rownames(self, i: str):
        self.__assayData = {key: x.loc[i] for key, x in self.__assayData.items()}
        self.__rowDat = self.rowData().loc[i]
        return self

    def __columns(self, i: str):      
        self.__assayData = {key: x[i] for key, x in self.__assayData.items()}
        self.__colDat = self.colData().loc[i]
        return self

    def setIndex(self, index):
        if len(index) == len(self.index()):
            self.__rowDat.index = index
            for x in self.__assayData.keys():
                self.__assayData[x].index = index

    def index(self):
        return self.rowData().index
    
    def setColnames(self, names):
        if len(names) == len(self.colnames()):
            self.__colDat.index = names
            for x in self.__assayData.keys():
                self.__assayData[x].columns = names
                

    
    def colnames(self):
        return self.colData().index

    def __setitem__(self, i, assay: pd.DataFrame):
        # check if assay dimnames == assayData dimnames
        self.assayData[i] = assay
    
    def __parseArgs(self, arg, lookup):
        if type(arg) == pd.core.series.Series:
            arg = np.where(arg)[0].tolist()
        if type(arg) == slice:
            arg = lookup[arg]
        elif type(arg) != list:
            if type(arg) == int:
                arg = lookup[arg]
            arg = [arg]
        elif type(arg) == list:
            for j in range(len(arg)):
                if type(arg[j]) == int:
                    arg[j] = lookup[arg[j]]
        return arg

    def __getitem__(self, i):
        try:
            rows, cols = i        
        except ValueError:
            rows = slice(None)
            cols = i
        
        rows = self.__parseArgs(rows, self.index())
        cols = self.__parseArgs(cols, self.colnames())
                
        
        if type(rows) != list and type(rows) != pd.core.indexes.base.Index:
            rows = [rows]
        if type(cols) != list and type(cols) != pd.core.indexes.base.Index:
            cols = [cols]    
        
        return copy.copy(self).__rownames(rows).__columns(cols)

    def assays(self):
        return self.__assayData

    def assay(self, i = 0):
        if type(i) == int:
            i = list(self.__assayData.keys())[i]
        return self.__assayData[i]
   
    def colData(self, col = None):
        if col == None:
            col = self.__colDat.columns
        return self.__colDat[col]
    
    def rowData(self, col = None):
        if col == None:
            col = self.__rowDat.columns
        return self.__rowDat[col]

    def metadata(self, key = None):
        if key == None:
            key = self.__metadat.keys()
        return self.__metadat[key]
    
    def __repr__(self) -> str:
        s = ["class: pySummarizedExperiment",
        "dim: %s %s" % (len(self.index()), len(self.colnames())),
        "metadata(%d): %s" % (len(self.__metadat), self.__metadat), 
        "rownames(%d): %s" % (len(self.index()), ", ".join(self.index()[0:5].astype(str))),
        "rowData names(%d): %s" % (len(self.rowData().columns),  ", ".join(self.rowData()[0:5].astype(str))),
        "colnames(%d): %s" % (len(self.colnames()), ", ".join(self.colnames()[0:5].astype(str))),
        "colData names(%d): %s" % (len(self.colData().columns), ", ".join(self.colData()[0:5].astype(str))),
        "assays(%d): %s" % (len(self.assays()), ", ".join(list(self.assays().keys())[0:5])) 
        ]
        return "\n".join(s)
    


