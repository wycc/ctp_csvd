
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

#### CTP DICOM 中的 InstanceNumber 說明

在 CT 灌注（CTP）影像中，**`InstanceNumber`** 是一個關鍵的 DICOM 標籤，代表：

1. **切片與時間的複合編號**：在 CTP 掃描中，`InstanceNumber` 通常按照以下模式編號：
   - 第一個時間點的所有切片：1, 2, 3, ..., N
   - 第二個時間點的所有切片：N+1, N+2, N+3, ..., 2N
   - 第三個時間點的所有切片：2N+1, 2N+2, 2N+3, ..., 3N
   - 以此類推...

2. **空間位置指示**：在同一時間點內，`InstanceNumber` 的順序通常對應切片的空間位置（從頭到腳或從腳到頭）

3. **時間序列標識**：透過 `InstanceNumber` 可以推算出：
   - 該切片屬於哪個時間點
   - 該切片在該時間點中的空間位置

**範例**：
假設 CTP 掃描有 20 個切片，5 個時間點：
- InstanceNumber 1-20：第 1 個時間點的切片
- InstanceNumber 21-40：第 2 個時間點的切片
- InstanceNumber 41-60：第 3 個時間點的切片
- InstanceNumber 61-80：第 4 個時間點的切片
- InstanceNumber 81-100：第 5 個時間點的切片

這就是為什麼 `readCTPSeries` 函數中使用 `fixedInterval` 參數時，會用以下公式計算時間點：
```python
df['refTime'] = np.floor((df.instanceNumber - 1).astype('float') / fixedInterval)
```
其中 `fixedInterval` 就是每個時間點的切片數量。

### 3. 時間參考設定階段 (第 74-91 行)

#### 3.1 使用實際時間戳記 (`refTime` 不為 None)

```python
if refTime is not None:
    if refTime == 'AcquisitionTime':
        df['refTime'] = df.acqTime
    elif refTime == 'ContentTime':
        df['refTime'] = df.contentTime
    else:
        raise NotImplementedError(f"{refTime} as refTime is not supported.")
    df['timeDiff'] = df.refTime.diff()
    df['timeDiff'].iloc[0] = pd.to_timedelta(0)
    locs = np.concatenate([[0], df[df.timeDiff > pd.to_timedelta('0.5s')].index])
```

- 根據指定的 `refTime` 設定參考時間
- 計算相鄰影像間的時間差
- 當時間差大於 0.5 秒時，認為是新的時間點，記錄其位置

#### 3.2 使用編號系統 (`refTime` 為 None)

```python
else:
    if fixedInterval is not None:
        df['refTime'] = np.floor((df.instanceNumber - 1).astype('float') / fixedInterval)
    else:
        df['refTime'] = df.acqNumber
    df['timeDiff'] = df.refTime.diff()
    df['timeDiff'].iloc[0] = 0
    locs = np.concatenate([[0], df[df.timeDiff > 0].index])
```

- 如果指定 `fixedInterval`，則根據實例編號計算時間點
- 否則使用 `AcquisitionNumber` 作為時間參考
- 當編號差異大於 0 時，認為是新的時間點

### 4. 影像重建階段 (第 93-109 行)

```python
times = []
images = []
for tp in range(len(locs)-1):
    st = locs[tp]
    ed = locs[tp+1]
    im = dcmLstToSitkImage(df.iloc[st:ed].dcmObj.tolist())
    times.append(df.iloc[st].refTime)
    images.append(im)
```

- 根據時間點分割位置 (`locs`) 將 DICOM 切片分組
- 每組切片使用 [`dcmLstToSitkImage`](ctp_csvd/ctp_tools/ctp.py:11) 函數重建為 3D 影像
- 記錄每個時間點的時間戳記

### 5. 時間正規化階段 (第 101-104 行)

```python
if refTime is not None:
    times = np.array([(i - times[0]).total_seconds() for i in times])
else:
    times = np.array([(i - times[0]) for i in times])
```

- 如果使用實際時間戳記，將時間轉換為相對於第一個時間點的秒數
- 如果使用編號系統，直接計算相對差值

## 關鍵特性

### 1. 自動時間點檢測
函數能夠自動檢測 CT 灌注序列中的不同時間點，無論是基於實際時間戳記還是編號系統。

### 2. 彈性的時間參考
支援多種時間參考方式：
- **AcquisitionTime**: 實際擷取時間
- **ContentTime**: 內容時間
- **AcquisitionNumber**: 擷取編號
- **固定間隔**: 基於實例編號的固定間隔分組

### 3. 強健的分組邏輯
- 使用時間差閾值（0.5秒）或編號差異來分組
- 確保每個時間點包含完整的 3D 體積

### 4. 除錯支援
提供除錯模式，返回詳細的中間處理結果，便於問題診斷。

## 使用場景

此函數主要用於：
1. **CT 灌注影像分析**: 處理動態 CT 掃描序列
2. **時間序列重建**: 將散亂的 DICOM 檔案重組為時間序列
3. **醫學影像預處理**: 為後續的灌注參數計算準備資料

## 相依性

函數依賴以下模組和函數：
- `pydicom`: DICOM 檔案讀取
- `pandas`: 資料結構化和處理
- `numpy`: 數值計算
- [`dcmLstToSitkImage`](ctp_csvd/ctp_tools/ctp.py:11): 將 DICOM 列表轉換為 SimpleITK 影像

## 注意事項

1. **InstanceNumber 的重要性**:
   - 函數假設 DICOM 檔案的 `InstanceNumber` 能正確反映切片的時空順序
   - 在 CTP 中，`InstanceNumber` 是切片位置和時間點的複合編號
   - 如果 `InstanceNumber` 不連續或不按預期順序，可能導致影像重建錯誤

2. **時間精度**: 時間差閾值設定為 0.5 秒，可能需要根據具體掃描協議調整

3. **記憶體使用**: 所有 DICOM 檔案會同時載入記憶體，大型資料集可能需要考慮記憶體限制

4. **固定間隔參數**: 當使用 `fixedInterval` 時，必須確保該值等於每個時間點的切片數量

## 範例使用

```python
# 基本使用
dcm_files = ['file1.dcm', 'file2.dcm', ...]
images, times = readCTPSeries(dcm_files)

# 使用 AcquisitionTime 作為時間參考
images, times = readCTPSeries(dcm_files, refTime='AcquisitionTime')

# 除錯模式
images, times, df, locs = readCTPSeries(dcm_files, debug=True)
```

這個函數是 CT 灌注影像處理流程中的關鍵組件，為後續的灌注參數計算（如 CBF、CBV、MTT 等）提供了標準化的輸入格式。
