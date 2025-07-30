## runTTP 函數分析


# TTP (Time to Peak) 到達峰值時間詳解

## 概述

TTP (Time to Peak) 是 **到達峰值時間** 的縮寫，是 CT 灌注成像 (CTP) 中最重要且最直觀的血流動力學參數之一。它反映了對比劑從注射到在特定組織中達到最高濃度所需的時間。

---

## 基本定義

### TTP 的含義
TTP 是指從對比劑開始注射到組織內對比劑濃度達到最高峰值所需的時間，通常以秒為單位測量。

### 物理意義
TTP 反映了對比劑在血管系統中的傳輸時間，包括：
- 從注射部位到心臟的時間
- 心臟泵血到主動脈的時間  
- 經由動脈系統到達目標組織的時間
- 在組織中達到平衡濃度的時間

---

## 血流傳輸路徑與時間

### 對比劑傳輸過程
```
注射部位 → 靜脈系統 → 右心 → 肺循環 → 左心 → 主動脈 → 頸動脈 → 腦動脈 → 組織
    ↓         ↓        ↓      ↓       ↓       ↓        ↓        ↓       ↓
  t=0秒     t=1秒    t=2秒   t=3秒   t=4秒   t=5秒    t=6秒    t=7秒   t=TTP
```

### 時間組成
```python
總 TTP = 循環時間 + 傳輸時間 + 組織平衡時間
       = (心臟+肺) + (動脈系統) + (毛細血管+組織)
       = 3-5秒    + 2-4秒     + 5-15秒
```

---

## 不同組織的 TTP 特徵

### 1. 動脈血管 (Arterial Vessels)

#### 時間特徵
- **TTP 範圍**: 4-8 秒
- **標準值**: 6 ± 2 秒
- **變異係數**: 低 (< 15%)

#### 生理基礎
- 直接接受動脈血流
- 血流速度最快
- 對比劑濃度變化最敏感

#### 臨床意義
```python
動脈 TTP 應用:
- 血管狹窄檢測: TTP 延長提示上游狹窄
- 側支循環評估: 多條動脈的 TTP 比較
- 心輸出量評估: 整體 TTP 延長提示心功能不全
```

### 2. 靜脈血管 (Venous Vessels)

#### 時間特徵
- **TTP 範圍**: 8-15 秒
- **標準值**: 12 ± 3 秒
- **變異係數**: 中等 (15-25%)

#### 生理基礎
- 接受組織回流血液
- 血流速度中等
- 對比劑濃度較動脈延遲

#### 臨床意義
```python
靜脈 TTP 應用:
- 靜脈回流評估: TTP 過度延長提示回流障礙
- 顱內壓評估: 靜脈 TTP 與顱內壓相關
- 血管畸形檢測: 動靜脈瘻會改變 TTP 模式
```

### 3. 腦組織 (Brain Tissue)

#### 時間特徵
- **TTP 範圍**: 15-25 秒
- **標準值**: 20 ± 5 秒
- **變異係數**: 高 (25-35%)

#### 生理基礎
- 經由毛細血管床灌注
- 血腦屏障影響對比劑分布
- 組織容積效應

#### 臨床意義
```python
腦組織 TTP 應用:
- 缺血檢測: TTP > 25秒 提示缺血
- 梗死診斷: TTP > 35秒 或無峰值提示梗死
- 灌注評估: TTP 分布圖顯示血流模式
```


### 函數簽名
```python
def runTTP(imCTC, tIdx, s0Idx, brainMask, outsideValue=-1):
```

### 功能描述
計算到達峰值時間 (Time to Peak, TTP)，這是評估血流動力學的重要參數。

### 參數說明
- **`imCTC`**: 對比劑時間曲線影像（來自 getCTC 函數的輸出）
- **`tIdx`**: 時間索引陣列
- **`s0Idx`**: 對比劑開始注射的時間點索引
- **`brainMask`**: 腦部遮罩
- **`outsideValue`**: 腦部外區域的填充值，預設為-1

### 主要處理流程

#### 1. 時間軸調整
首先，我們算的時候只考慮CTC 偵測到顯影劑後開始算起。s0Idx 之前的影像都不算。

```python
tIdxFil = tIdx[s0Idx:] - tIdx[s0Idx]  # 重設時間起點為0
imCTCFil = imCTC[s0Idx:]              # 對應的影像序列
```
- 以對比劑注射時間為起點重新計算時間軸
- 提取注射後的影像序列

#### 2. 峰值時間計算
接未來就每一個點用 argmax 找出最大值的是第幾個 3D 影像。

```python
TTP = tIdxFil[imCTCFil.argmax(axis=0)]
```
- 找到每個像素在時間軸上的最大值位置
- 對應的時間即為該像素的到達峰值時間

TTP 中存的就是每一點最大值的時間位置。

#### 3. 遮罩處理
我們把不是腦的部份移除。如果去頭骨己經在前處理做了，這一步就是不必要的。

```python
TTP[sitk.GetArrayFromImage(brainMask) == 0] = outsideValue
```
- 將腦部外的區域設為指定的填充值
- 確保只有腦部區域有有效的 TTP 值

### 返回值
- **`TTP`**: 每個像素的到達峰值時間圖

---

## 函數間的關聯性

### 工作流程
1. **getCTC** 函數負責前處理：
   - 從原始影像序列中提取對比劑濃度變化
   - 檢測對比劑注射起始時間
   - 進行基線校正

2. **runTTP** 函數使用 getCTC 的輸出：
   - 計算每個像素的對比劑到達峰值時間
   - 生成 TTP 參數圖

### 臨床意義
- **TTP** 是重要的血流動力學參數
- 在腦中風診斷中，TTP 延遲可能表示血流灌注不足
- 常與其他參數（CBF, CBV, MTT）一起用於評估腦組織的血流狀態

---

## 技術特點

### 優點
1. **自動化檢測**: 自動檢測對比劑注射時間點
2. **基線校正**: 有效去除基線影響
3. **遮罩處理**: 確保只分析腦部區域
4. **參數化設計**: 可調整閾值參數適應不同情況

### 注意事項
1. **閾值敏感性**: HU 閾值和比例閾值需要根據具體情況調整
2. **時間解析度**: 需要足夠的時間解析度來準確檢測峰值
3. **影像品質**: 影像雜訊可能影響峰值檢測的準確性

---

---

## runAIF 函數分析
一般而言，我們可以手動標出主要動脈的位置。把它的亮度當成是一個期準值來用。這個函數基本上是自動找出來做適合當基準值的位置。這樣一來
我們就不用一個一個標注主動脈的位置了。它的原理是用下面三個條件來決定
- 基於三個條件篩選候選區域：
  - 曲線下面積 (AUC) 大於閾值
  - TTP 小於閾值（早期到達峰值）
  - TTP 大於0（有效區域）

基本上就是在最早 TTP 達到峰值的位置中，使用形態學的方法找一個符合所有條件的位置來用。詳細的條件請看下面分析。

### 函數簽名
```python
def runAIF(imCTC, tIdx, brainMask, TTP,
           roi=None,
           TTPThres=8, AUCThres=1000,
           dilateRadius=5, erodeRadius=2,
           candidateVolThres = 2000
          ):
```

### 功能描述
自動識別動脈輸入函數 (Arterial Input Function, AIF)，這是 CTP 分析中的關鍵步驟。AIF 代表對比劑在動脈中的濃度時間曲線，用於後續的去卷積分析。

### 參數說明
- **`imCTC`**: 對比劑時間曲線影像
- **`tIdx`**: 時間索引陣列
- **`brainMask`**: 腦部遮罩
- **`TTP`**: 到達峰值時間圖
- **`roi`**: 感興趣區域 (可選)，格式為 [x_min, y_min, z_min, x_max, y_max, z_max]
- **`TTPThres`**: TTP 閾值，預設為8秒
- **`AUCThres`**: 曲線下面積閾值，預設為1000
- **`dilateRadius`**: 膨脹半徑，預設為5像素
- **`erodeRadius`**: 侵蝕半徑，預設為2像素
- **`candidateVolThres`**: 候選區域體積閾值，預設為2000立方毫米

### 主要處理流程

#### 1. 初始候選區域識別
```python
brainMaskNp = sitk.GetArrayFromImage(brainMask)
aifCand = (imCTC.sum(0) * brainMaskNp > AUCThres) * (TTP < TTPThres) * (TTP > 0)
```
- 基於三個條件篩選候選區域：
  - 曲線下面積 (AUC) 大於閾值
  - TTP 小於閾值（早期到達峰值）
  - TTP 大於0（有效區域）

#### 2. ROI 限制（可選）
```python
if roi is not None:
    aifCandRoi = np.zeros_like(aifCand)
    aifCandRoi[roi[2]:roi[5], roi[1]:roi[4], roi[0]:roi[3]] = aifCand[roi[2]:roi[5], roi[1]:roi[4], roi[0]:roi[3]]
    aifCand = aifCandRoi
```
- 如果指定了 ROI，則將候選區域限制在該範圍內

#### 3. 形態學處理
```python
aifCand = binary_erosion(
    binary_dilation(aifCand, np.ones([1, dilateRadius, dilateRadius], bool)),
    np.ones([1, erodeRadius, erodeRadius], bool))
```
- 先膨脹後侵蝕（閉運算）
- 去除小的雜訊區域，連接鄰近的候選區域

#### 4. 連通組件分析
```python
aifCand, nComp = cc3d.connected_components(aifCand, return_N=True)
```
- 識別獨立的連通組件
- 每個組件代表一個潛在的 AIF 候選區域

#### 5. 候選區域評估
```python
for idx in range(1, nComp+1):
    curCand = aifCand == idx
    vol = curCand.sum()
    curve = imCTC[:, curCand]
    curveMean = curve.mean(axis=1)
    peakEndDiff = np.max(curveMean) - np.mean(curveMean[-3:])
```

對每個候選區域進行評估：
- **體積檢查**: 排除過大的區域（可能是靜脈）
- **曲線擬合**: 使用 gamma variate 函數擬合時間曲線
- **誤差計算**: 評估擬合品質

#### 6. Gamma Variate 函數擬合
```python
try:
    popts, _ = fit_gv(tIdx, curveMean)
except:
    traceback.print_exc()
    continue
err = np.sqrt(np.sum((curve - gv(tIdx, *popts)[:, np.newaxis]) ** 2, axis=1))
```
- 使用 gamma variate 函數擬合平均時間曲線
- 計算擬合誤差

#### 7. 候選區域評分
```python
cands = pd.DataFrame(cands)
cands['score'] = cands.vol * cands.peakEndDiff / cands.meanErr
bestCand = np.argmax(cands.score)
```
- 綜合評分公式：`體積 × 峰值差異 / 平均誤差`
- 選擇評分最高的候選區域作為最佳 AIF

### 返回值
- **`aifProps`**: 最佳 AIF 的 gamma variate 參數
- **`aifCand==aifSegIdx`**: 最佳 AIF 區域的二值遮罩
- **`cands`**: 所有候選區域的詳細資訊 DataFrame

### 評分機制詳解

#### 評分因子
1. **體積 (vol)**:
   - 較大的血管通常有更好的信噪比
   - 但過大可能是靜脈而非動脈

2. **峰值差異 (peakEndDiff)**:
   - `np.max(curveMean) - np.mean(curveMean[-3:])`
   - 動脈應該有明顯的峰值和快速的洗出

3. **平均誤差 (meanErr)**:
   - 擬合品質的指標
   - 較小的誤差表示更好的 gamma variate 擬合

#### 最終評分
```python
score = vol * peakEndDiff / meanErr
```
- 偏好體積適中、峰值明顯、擬合良好的區域
- 這個評分機制有助於區分動脈和靜脈

### 技術特點

#### 優點
1. **自動化識別**: 無需手動選擇 AIF 區域
2. **多重篩選**: 結合 AUC、TTP、體積等多個條件
3. **形態學處理**: 有效去除雜訊和小區域
4. **定量評分**: 客觀的評分機制選擇最佳候選

#### 限制
1. **參數敏感性**: 多個閾值參數需要根據具體情況調整
2. **解剖變異**: 可能無法處理異常的血管解剖結構
3. **影像品質依賴**: 需要足夠的時間和空間解析度

### 臨床意義
- **AIF 品質**: 直接影響後續 CBF、CBV、MTT 等參數的準確性
- **動脈選擇**: 通常選擇大腦中動脈或前動脈作為 AIF
- **病理影響**: 血管狹窄或閉塞可能影響 AIF 的自動識別

---

---

