# Example
Start by importing the libraries


```python
import pandas as pd
from pySummarizedExperiment import pySummarizedExperiment
```

## Using predefined Pandas DataFrames
Using three dataframes, we can create a pySummarizedExperiment. 


```python
data = [[1,2,3,4], [5,6,7,8], [6, 7, 8, 9]]
columns = ["a", "b", "c", "d"]
rownames = ["FT1", "FT2", "FT3"]
data = pd.DataFrame(data, columns = columns, index = rownames)
rowData = pd.DataFrame([1,2,3], index = rownames, columns = ["mz"])
colData = pd.DataFrame([["QC", 1], ["QC", 2], ["BLANK", 3], ["BLANK", 4]],  index=columns, columns = ["Type", "Injection"])
assays = {"first_assay": data, "second_assay": data * 2}
exp = pySummarizedExperiment(assays = assays, colData = colData, rowData = rowData)
exp
```




    class: pySummarizedExperiment
    dim: 3 4
    metadata(1): Datetime
    rownames(3): FT1, FT2, FT3
    rowData names(1): mz
    colnames(4): a, b, c, d
    colData names(2): Type, Injection
    assays(2): first_assay, second_assay



## From a long-format DataFrame
We can also create a pySummarizedExperiment using a single long-format dataframe. We need a rowIndex and colIndex to define how to create the rowData and Coldata. Next, the cardinality between columns is checked to define how columns should be assigned. First we create a long dataframe.


```python
from random import randint
from random import seed

seed(42)
df = pd.DataFrame({
    "samples": ["a", "b", "c", "d"] * 5,
    "sample_day": [1, 4, 4, 7] * 5,
    "features": [1, 2, 3, 4, 5] * 4,
    "feature_polarity": ["pos", "neg", "neg", "pos", "pos"] * 4,
    "assay": [randint(0, 1000) for _ in range(20)],
    "assay_other": [_ for _ in range(20)]
})
df.head(5)

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>samples</th>
      <th>sample_day</th>
      <th>features</th>
      <th>feature_polarity</th>
      <th>assay</th>
      <th>assay_other</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>a</td>
      <td>1</td>
      <td>1</td>
      <td>pos</td>
      <td>654</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>b</td>
      <td>4</td>
      <td>2</td>
      <td>neg</td>
      <td>114</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>c</td>
      <td>4</td>
      <td>3</td>
      <td>neg</td>
      <td>25</td>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>d</td>
      <td>7</td>
      <td>4</td>
      <td>pos</td>
      <td>759</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>a</td>
      <td>1</td>
      <td>5</td>
      <td>pos</td>
      <td>281</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>



Next, we create a `pySummarizedExperiment` by setting the dataframe as the `longDf` parameter. Here we use the features column as `rowIndex` parameter and the samples column as `colIndex`. When we print the object, we can see the following:
* The class name
* The dimensions (rows x columns) of the assay. 
* Rownames (row index) of the rowData and Assays
* rowData names, the columns of the rowData
* colnames, the columns of the assay and rowIndex of the colData
* colData names, columns of the colData
* assays, names of the stored assays, first one is used by default 



```python
exp = pySummarizedExperiment(longDf = df, rowIndex = "features", colIndex = "samples")
print(exp)
```

    class: pySummarizedExperiment
    dim: 5 4
    metadata(1): Datetime
    rownames(5): 1, 2, 3, 4, 5
    rowData names(1): feature_polarity
    colnames(4): a, b, c, d
    colData names(1): sample_day
    assays(2): assay, assay_other
    

## Accessing Data
pySummarizedExperiment implements the basic functionality of colData, rowData, and assay functions from R as methods of the object. Here we showcase the methods `.coldata()`, `.rowData()` and `.assay()`. For colData and rowData, it is possible to give one or more column names to receive a subset. For assays, either an assay name or index number of the assay can be given te receive a specific assay. By default, the assay with the first index will be returned.


```python
# colData
exp.colData()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sample_day</th>
    </tr>
    <tr>
      <th>samples</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>a</th>
      <td>1</td>
    </tr>
    <tr>
      <th>b</th>
      <td>4</td>
    </tr>
    <tr>
      <th>c</th>
      <td>4</td>
    </tr>
    <tr>
      <th>d</th>
      <td>7</td>
    </tr>
  </tbody>
</table>
</div>




```python
# rowData
exp.rowData()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>feature_polarity</th>
    </tr>
    <tr>
      <th>features</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>pos</td>
    </tr>
    <tr>
      <th>2</th>
      <td>neg</td>
    </tr>
    <tr>
      <th>3</th>
      <td>neg</td>
    </tr>
    <tr>
      <th>4</th>
      <td>pos</td>
    </tr>
    <tr>
      <th>5</th>
      <td>pos</td>
    </tr>
  </tbody>
</table>
</div>




```python
# assay, either use an index or name, both are equal
exp.assay(0) 
exp.assay("assay")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>samples</th>
      <th>a</th>
      <th>b</th>
      <th>c</th>
      <th>d</th>
    </tr>
    <tr>
      <th>features</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>654</td>
      <td>250</td>
      <td>692</td>
      <td>604</td>
    </tr>
    <tr>
      <th>2</th>
      <td>432</td>
      <td>114</td>
      <td>228</td>
      <td>758</td>
    </tr>
    <tr>
      <th>3</th>
      <td>913</td>
      <td>32</td>
      <td>25</td>
      <td>142</td>
    </tr>
    <tr>
      <th>4</th>
      <td>754</td>
      <td>558</td>
      <td>30</td>
      <td>759</td>
    </tr>
    <tr>
      <th>5</th>
      <td>281</td>
      <td>104</td>
      <td>89</td>
      <td>95</td>
    </tr>
  </tbody>
</table>
</div>



## Subsetting
The `pySummarizedExperiment` object aims to have similar syntax as Pandas DataFrames. Therefore there are three ways of subsetting data. Using `lists`, `slices` or `masks`. Each of these will be covered below.
### Using lists
We can directly subset the object by providing lists with the indexes and columns. Here we are selecting the 2nd and 3rd row using indexes (1 and 2) and columns 'b' and 'c'. The result is a pySummarizedExperiment with two-by-two dimensions 


```python
exp[[1,2], ["b", "c"]]
```




    class: pySummarizedExperiment
    dim: 2 2
    metadata(1): Datetime
    rownames(2): 2, 3
    rowData names(1): feature_polarity
    colnames(2): b, c
    colData names(1): sample_day
    assays(2): assay, assay_other



### Using slices
We can subset the pySummarizedExperiment similar to Pandas syntax. Here we take all rows by supplying an empty slice `:` and take the samples `a` and `b` by supplying a slice of 0:2, which includes columns 0 and 1 (which is equal to columns 'a' and 'b'). Finally, we use the `.assay()` method to receive the assay with the subsetted object.


```python
exp[:, 0:2].assay()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>samples</th>
      <th>a</th>
      <th>b</th>
    </tr>
    <tr>
      <th>features</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>654</td>
      <td>250</td>
    </tr>
    <tr>
      <th>2</th>
      <td>432</td>
      <td>114</td>
    </tr>
    <tr>
      <th>3</th>
      <td>913</td>
      <td>32</td>
    </tr>
    <tr>
      <th>4</th>
      <td>754</td>
      <td>558</td>
    </tr>
    <tr>
      <th>5</th>
      <td>281</td>
      <td>104</td>
    </tr>
  </tbody>
</table>
</div>



### Using masks
A pySummarizedExperiment supports Pandas-style boolean masking for subsetting. Here we are subsetting rows by filtering based on the `feature_polarity` column in the rowData. We select only rows that have `pos` as value. We can see that we only retain rows 1, 4, and 5


```python
exp[exp.rowData("feature_polarity") == "pos", :]
```




    class: pySummarizedExperiment
    dim: 3 4
    metadata(1): Datetime
    rownames(3): 1, 4, 5
    rowData names(1): feature_polarity
    colnames(4): a, b, c, d
    colData names(1): sample_day
    assays(2): assay, assay_other



## Plotting
Since `pySummarizedExperiment` is based on Pandas DataFrames, it is easy to inspect data via plots. Some plots can be made directly using plots built in Pandas itself, whereas others can be used using other libraries as `matplotlib` or `seaborn`. Here, we construct a boxplot using the `.boxplot()` method. It accepts the following parameters:
* `assay` Which assay to plot
* `log2`: Transform the given assay with a log2 transformation
* `autoscale`: Transform the assay by scaling and centering
* `transpose`: If true, plot the rows as boxes, plot the columns otherwise (default)

The result can be seen underneath.


```python
exp.boxplot(assay = 0, log2 = True, autoscale=False, transpose = False)
```




    <AxesSubplot:>




    
![png](example_files/example_19_1.png)
    



```python
exp.boxplot(assay = 0, log2 = True, autoscale=True, transpose = True)
```




    <AxesSubplot:>




    
![png](example_files/example_20_1.png)
    


## Principle Component Analysis (PCA)
One of the most common transformations in -omics fields is a Principle Component Analysis (PCA). Similar to the `.boxplot()` method, this method accepts the `log2`, `autoscale` and `transpose` parameters. However, it also accepts the parameter `components`, which indicates the amount of components to create (defaults to 5). Unlike the `.boxplot()` method, this method returns the PCA DataFrame. Thus it can be plotted by using the `.plot.scatter()` in-built Pandas method.


```python
exp.pca(assay = "assay", components = 2, log2 = True, autoscale=True, transpose=True).plot.scatter(x = 0, y = 1)
```




    <AxesSubplot:xlabel='0', ylabel='1'>




    
![png](example_files/example_22_1.png)
    


## Export to long format
To export the data to a single dataframe, we can use the method `toLongDataFrame()`. It will return all data with the exception for metadata. This method is useful for writing data to a single `.csv` or similar data file. 


```python
exp.toLongDataFrame().head(5)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>features</th>
      <th>samples</th>
      <th>assay</th>
      <th>assay_other</th>
      <th>feature_polarity</th>
      <th>sample_day</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>a</td>
      <td>654</td>
      <td>0</td>
      <td>pos</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>a</td>
      <td>432</td>
      <td>16</td>
      <td>neg</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>a</td>
      <td>913</td>
      <td>12</td>
      <td>neg</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>a</td>
      <td>754</td>
      <td>8</td>
      <td>pos</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>a</td>
      <td>281</td>
      <td>4</td>
      <td>pos</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>


