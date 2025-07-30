# getCTC 和 runTTP 函數分析報告

## 概述

`getCTC` 和 `runTTP`。這兩個函數是 CTP (CT Perfusion) 灌注分析流程中的重要組件，用於處理醫學影像並計算血流動力學參數。

# CTC (Contrast Time Curve) 對比劑時間曲線詳解

## 概述

CTC (Contrast Time Curve) 是 **對比劑時間曲線** 的縮寫，是 CT 灌注成像 (CTP) 分析的核心概念。它描述了對比劑在組織中濃度隨時間變化的曲線，是評估血流動力學的重要工具。

---

## 基本概念

### 定義
對比劑時間曲線是指在 CT 灌注掃描過程中，特定組織或血管內對比劑濃度隨時間變化的函數關係。

### 物理基礎
- **對比劑**：通常使用含碘的造影劑（如碘海醇、碘普羅胺等）
- **注射方式**：經靜脈快速注射，通常使用高壓注射器
- **影像表現**：對比劑增加組織的 X 射線吸收係數，在 CT 影像中表現為 HU 值增加

---

## CTC 的典型形態

### 視覺化表示
```
HU 值 (亨氏單位)
 ↑
 │        ╭─╮ ← 峰值期
 │       ╱   ╲
 │      ╱     ╲
 │     ╱       ╲ ← 洗出期
 │    ╱         ╲
 │   ╱           ╲___
 │  ╱                ╲___
 │ ╱ ← 上升期              ╲___
 │╱                           ╲___
 └─────────────────────────────────→ 時間 (秒)
基線期  注射開始  到達峰值    洗出完成
```

### 四個主要階段

#### 1. 基線期 (Baseline Phase)
- **時間**：對比劑注射前
- **特徵**：HU 值保持穩定的低值
- **持續時間**：通常 5-10 秒
- **臨床意義**：建立組織的基礎密度值

#### 2. 上升期 (Upslope Phase)
- **時間**：對比劑開始到達組織
- **特徵**：HU 值快速上升
- **斜率**：反映血流速度和血管通透性
- **臨床意義**：評估組織灌注狀態

#### 3. 峰值期 (Peak Phase)
- **時間**：對比劑濃度達到最高點
- **特徵**：HU 值達到最大值
- **參數**：TTP (Time to Peak) - 到達峰值時間
- **臨床意義**：反映血流到達時間

#### 4. 洗出期 (Washout Phase)
- **時間**：對比劑逐漸被清除
- **特徵**：HU 值緩慢下降
- **斜率**：反映血流清除能力
- **臨床意義**：評估血管通透性和組織代謝

---

## 不同組織的 CTC 特徵

### 動脈血管 (Arterial Vessels)

#### 特徵
- **到達時間**：最早（通常 8-12 秒）
- **TTP**：短（< 8 秒）
- **峰值**：最高
- **上升斜率**：最陡峭
- **洗出速度**：快速

#### 典型參數
```python
動脈 CTC 參數：
- TTP: 4-8 秒
- 峰值 HU: 200-400
- 上升斜率: 高
- 洗出半衰期: 短
```

### 靜脈血管 (Venous Vessels)

#### 特徵
- **到達時間**：中等（通常 12-20 秒）
- **TTP**：中等（8-15 秒）
- **峰值**：中等
- **上升斜率**：中等
- **洗出速度**：中等

#### 典型參數
```python
靜脈 CTC 參數：
- TTP: 8-15 秒
- 峰值 HU: 100-200
- 上升斜率: 中等
- 洗出半衰期: 中等
```

### 腦組織 (Brain Tissue)

#### 特徵
- **到達時間**：最晚（通常 15-25 秒）
- **TTP**：長（> 15 秒）
- **峰值**：最低
- **上升斜率**：最緩慢
- **洗出速度**：最慢

#### 典型參數
```python
腦組織 CTC 參數：
- TTP: 15-25 秒
- 峰值 HU: 20-80
- 上升斜率: 低
- 洗出半衰期: 長
```

---

## CTC 的數學模型

### Gamma Variate 函數

最常用的 CTC 擬合模型是 Gamma Variate 函數：

```python
def gamma_variate(t, t0, alpha, beta):
    """
    Gamma Variate 函數
    
    參數:
    t: 時間
    t0: 對比劑到達時間 (秒)
    alpha: 上升斜率參數
    beta: 洗出速率參數 (秒)
    """
    return np.maximum(0, t-t0)**alpha * np.exp(-t/beta) * (np.sign(t - t0) + 1) / 2
```

#### 參數意義
- **t0 (到達時間)**：對比劑首次到達組織的時間
- **alpha (形狀參數)**：控制上升期的陡峭程度
- **beta (尺度參數)**：控制洗出期的速度

### 其他數學模型

#### 1. 指數衰減模型
```python
def exponential_decay(t, A, k, t0):
    return A * np.exp(-k * (t - t0)) * (t >= t0)
```

#### 2. 雙指數模型
```python
def double_exponential(t, A1, k1, A2, k2, t0):
    return (A1 * np.exp(-k1 * (t - t0)) + A2 * np.exp(-k2 * (t - t0))) * (t >= t0)
```

---

## CTC 分析在程式碼中的實現

### getCTC 函數解析

從提供的程式碼中，`getCTC` 函數執行以下步驟：

#### 1. 對比劑總量計算
```python
totContrastAgentVal = np.array([
    imNp[t][((imNp[t] > huThres) & brainMaskNp).astype('bool')].sum() 
    for t in range(imNp.shape[0])
]) + 1
```
- 計算每個時間點腦部區域內超過 HU 閾值的像素總和
- `huThres=150`：HU 值閾值，用於識別對比劑

#### 2. 基線標準化
```python
totContrastAgentVal = (totContrastAgentVal - totContrastAgentVal[0]) / totContrastAgentVal[0]
```
- 以第一個時間點為基線
- 計算相對於基線的變化比例

#### 3. 注射時間點檢測
```python
cands = np.where(totContrastAgentVal > ratioThres)[0]
s0Idx = cands[0]  # 對比劑開始注射的時間點
```
- `ratioThres=0.05`：比例閾值
- 自動檢測對比劑注射開始時間

#### 4. 基線減法
```python
s0ImNp = imNp[:s0Idx].mean(axis=0)  # 注射前影像平均
imNp = imNp - s0ImNp  # 減去基線
imNp[imNp<0] = 0  # 負值設為0
```
- 使用注射前影像作為基線
- 進行基線校正

---

## 血流動力學參數計算

### 基於 CTC 的參數

#### 1. TTP (Time to Peak)
```python
def calculate_TTP(ctc, time_points):
    return time_points[np.argmax(ctc)]
```

#### 2. Peak Enhancement
```python
def calculate_peak(ctc):
    return np.max(ctc)
```

#### 3. AUC (Area Under Curve)
```python
def calculate_AUC(ctc, time_points):
    return np.trapz(ctc, time_points)
```

#### 4. Upslope
```python
def calculate_upslope(ctc, time_points):
    peak_idx = np.argmax(ctc)
    return ctc[peak_idx] / time_points[peak_idx]
```

### 進階參數（需要去卷積分析）

#### 1. CBF (Cerebral Blood Flow)
- **定義**：腦血流量
- **單位**：ml/100g/min
- **計算**：基於殘留函數的最大值

#### 2. CBV (Cerebral Blood Volume)
- **定義**：腦血容量
- **單位**：ml/100g
- **計算**：基於殘留函數的積分

#### 3. MTT (Mean Transit Time)
- **定義**：平均通過時間
- **單位**：秒
- **計算**：CBV/CBF

---

## 臨床應用

### 急性中風診斷

#### 正常組織 vs 缺血組織

| 參數 | 正常組織 | 缺血組織 | 梗死組織 |
|------|----------|----------|----------|
| TTP | 正常 | 延長 | 顯著延長或無峰值 |
| 峰值 | 正常 | 降低 | 顯著降低或無 |
| CBF | 正常 | 降低 | 顯著降低 |
| CBV | 正常 | 可能增加（側支循環） | 顯著降低 |

#### 缺血半暗帶識別
```python
# 缺血半暗帶的典型特徵
penumbra_criteria = {
    'TTP': '延長 > 4秒',
    'CBF': '降低 < 30% 正常值',
    'CBV': '相對保持或輕度降低',
    'MTT': '顯著延長'
}
```

### 腫瘤評估

#### 血管新生評估
- **高血管化腫瘤**：早期強化，高峰值
- **低血管化腫瘤**：延遲強化，低峰值
- **壞死區域**：無明顯強化

#### 治療反應監測
- **有效治療**：血管化減少，CTC 峰值降低
- **無效治療**：血管化維持或增加

---

## 技術考量

### 影像獲取參數

#### 時間解析度
- **推薦**：1-2 秒間隔
- **最小**：0.5 秒間隔
- **影響**：時間解析度影響 TTP 測量精度

#### 空間解析度
- **推薦**：512×512 矩陣
- **層厚**：5-10 mm
- **影響**：空間解析度影響小血管檢測

#### 掃描持續時間
- **推薦**：60-90 秒
- **最小**：45 秒
- **影響**：需要涵蓋完整的洗出期

### 對比劑注射參數

#### 注射速度
- **推薦**：4-6 ml/s
- **影響**：注射速度影響峰值和上升斜率

#### 對比劑劑量
- **推薦**：1-1.5 ml/kg
- **影響**：劑量影響信噪比和峰值

#### 注射延遲
- **推薦**：5-8 秒延遲後開始掃描
- **影響**：確保捕捉到基線期

---

## 品質控制

### 常見問題和解決方案

#### 1. 運動偽影
- **問題**：患者移動導致 CTC 異常
- **解決**：影像配準、固定裝置

#### 2. 對比劑外滲
- **問題**：注射部位對比劑滲漏
- **解決**：檢查注射部位、重新注射

#### 3. 心律不整
- **問題**：心律不整影響對比劑到達
- **解決**：心電圖監測、延長掃描時間

#### 4. 血管狹窄
- **問題**：上游血管狹窄影響 CTC 形態
- **解決**：結合血管造影、調整分析參數

### 品質評估指標

#### 1. 信噪比 (SNR)
```python
def calculate_SNR(ctc, baseline_std):
    peak_value = np.max(ctc)
    return peak_value / baseline_std
```

#### 2. 對比雜訊比 (CNR)
```python
def calculate_CNR(ctc_tissue, ctc_background, noise_std):
    return abs(np.max(ctc_tissue) - np.max(ctc_background)) / noise_std
```

#### 3. 擬合品質 (R²)
```python
def calculate_r_squared(observed, fitted):
    ss_res = np.sum((observed - fitted) ** 2)
    ss_tot = np.sum((observed - np.mean(observed)) ** 2)
    return 1 - (ss_res / ss_tot)
```

---

## 未來發展方向

### 技術改進

#### 1. 人工智慧應用
- **深度學習**：自動 CTC 分析和參數提取
- **機器學習**：異常 CTC 模式識別
- **神經網路**：去雜訊和影像增強

#### 2. 新的數學模型
- **分數階微積分**：更精確的 CTC 建模
- **多室模型**：考慮組織異質性
- **貝葉斯方法**：不確定性量化

#### 3. 多模態整合
- **CTP + MRP**：結合 CT 和 MR 灌注
- **CTP + DSA**：結合血管造影
- **CTP + PET**：結合代謝資訊

### 臨床應用擴展

#### 1. 新的疾病領域
- **心肌灌注**：冠心病評估
- **腎臟灌注**：腎功能評估
- **肝臟灌注**：肝病診斷

#### 2. 治療指導
- **個人化治療**：基於 CTC 的治療選擇
- **療效預測**：CTC 參數預測治療反應
- **預後評估**：長期預後預測

---

## 總結

CTC (Contrast Time Curve) 是 CT 灌注成像的核心概念，通過分析對比劑在組織中的時間變化模式，可以獲得重要的血流動力學資訊。理解 CTC 的基本原理、數學模型和臨床應用，對於正確解釋 CTP 結果和指導臨床決策至關重要。

隨著技術的不斷發展，CTC 分析將變得更加精確和自動化，為臨床診斷和治療提供更有價值的資訊。

---

## getCTC 函數分析

### 函數簽名
```python
def getCTC(imLst, brainMask, huThres=150, ratioThres=0.05):
```

### 功能描述
計算對比劑時間曲線 (Contrast Time Curve, CTC)，這是 CTP 分析的核心前處理步驟。

### 參數說明
- **`imLst`**: 影像序列列表 (SimpleITK Image 物件列表)
- **`brainMask`**: 腦部遮罩 (SimpleITK Image 物件)
- **`huThres`**: HU值閾值，預設為150，用於識別對比劑
- **`ratioThres`**: 比例閾值，預設為0.05，用於檢測對比劑注射時間點

### 主要處理流程

#### 1. 影像資料轉換
```python
imNp = np.stack([sitk.GetArrayFromImage(im) for im in imLst])
brainMaskNp = sitk.GetArrayFromImage(brainMask)
```
- 將 SimpleITK 影像列表轉換為 numpy 陣列
- 轉換腦部遮罩為 numpy 陣列

#### 2. 對比劑總量計算
第一步是要找出顯影劑打入的時間點。我們可以用整體亮度超過 huThres 的點數來決定。還沒有打顯影劑前都會在一個低的水平。

先把每一組 3D 影像中超過 threshold 的點數算出來。

```python
totContrastAgentVal = np.array(
    [imNp[t][((imNp[t] > huThres) & brainMaskNp).astype('bool')].sum() for t in range(imNp.shape[0])]
) + 1
```
- 計算每個時間點腦部區域內超過 HU 閾值的像素總和
- 加1避免除零錯誤

#### 3. 基線標準化
把第一組 3D 影像當成是基準點。把它減掉後再除以這個基準值。這樣可以轉換成一個相對的亮度。

```python
totContrastAgentVal = (totContrastAgentVal - totContrastAgentVal[0]) / totContrastAgentVal[0]
```
- 以第一個時間點為基線
- 計算相對於基線的變化比例

#### 4. 對比劑注射時間點檢測
```python
cands = np.where(totContrastAgentVal > ratioThres)[0]
if cands.shape[0] > 0:
    s0Idx = cands[0]
```
- 找到對比劑濃度超過閾值的第一個時間點
- 此時間點被認定為對比劑開始注射的時間

#### 5. 基線減法處理
```python
s0ImNp = imNp[:s0Idx].mean(axis=0)  # 注射前影像平均
imNp = imNp - s0ImNp                # 減去基線
imNp[imNp<0] = 0                    # 負值設為0
```
- 計算注射前影像的平均值作為基線
- 從所有影像中減去基線
- 將負值設為0，確保對比劑濃度為非負值

### 返回值
- **`imNp`**: 處理後的對比劑濃度影像序列
- **`s0Idx`**: 對比劑開始注射的時間點索引

### 錯誤處理
如果未檢測到對比劑注射（即沒有時間點超過閾值），函數返回 `(None, None)`。

---

