
# readCTPSeries 函數分析說明

## 概述

[`readCTPSeries`](ctp_csvd/ctp_tools/ctp.py:58) 是一個用於讀取和處理 CT 灌注（CTP, CT Perfusion）影像序列的函數。此函數能夠將 DICOM 檔案列表轉換為時間序列的 3D 影像，並根據時間資訊進行適當的分組和排序。

## 函數簽名

```python
def readCTPSeries(dcmList, refTime=None, fixedInterval=None, debug=False):
```

## 參數說明

- **`dcmList`**: DICOM 檔案路徑的列表
- **`refTime`**: 參考時間類型，可選值：
  - `None`: 使用 AcquisitionNumber 作為時間參考
  - `'AcquisitionTime'`: 使用 AcquisitionTime 作為時間參考
  - `'ContentTime'`: 使用 ContentTime 作為時間參考
- **`fixedInterval`**: 固定間隔參數，當 `refTime=None` 時使用
- **`debug`**: 是否返回除錯資訊

## 返回值

- **正常模式** (`debug=False`): 
  - `images`: 3D 影像列表
  - `times`: 對應的時間陣列
- **除錯模式** (`debug=True`):
  - `images`: 3D 影像列表
  - `times`: 對應的時間陣列
  - `df`: 包含所有 DICOM 資訊的 DataFrame
  - `locs`: 時間點分割位置

## 工作原理詳解

### 1. DICOM 檔案讀取階段 (第 59-61 行)

```python
ds = []
for i in dcmList:
    ds.append(pydicom.read_file(i))
```

函數首先使用 `pydicom.read_file()` 讀取所有 DICOM 檔案，將其儲存在 `ds` 列表中。

### 2. 資料結構化階段 (第 62-73 行)

```python
df = []
for i in ds:
    df.append({
        'acqTime': pd.to_datetime(f'{i.AcquisitionDate} {i.AcquisitionTime}'),
        'contentTime': pd.to_datetime(f'{i.ContentDate} {i.ContentTime}'),            
        'instanceNumber': i.InstanceNumber,
        'acqNumber': i.AcquisitionNumber,
        'dcmObj': i
    })
df = pd.DataFrame(df)
df.sort_values(by='instanceNumber', inplace=True)
df.reset_index(drop=True, inplace=True)
```

此階段將 DICOM 資訊轉換為結構化的 pandas DataFrame，包含：
- **`acqTime`**: 擷取時間（AcquisitionDate + AcquisitionTime）
- **`contentTime`**: 內容時間（ContentDate + ContentTime）
- **`instanceNumber`**: 實例編號（詳見下方說明）
- **`acqNumber`**: 擷取編號
- **`dcmObj`**: 原始 DICOM 物件

資料按 `instanceNumber` 排序以確保正確的順序。
