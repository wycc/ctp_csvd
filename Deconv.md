## runDeconv 函數分析
#### 去卷積原理
在 CTP 分析中，組織濃度曲線 C(t) 可以表示為：
```
C(t) = CBF × AIF(t) ⊗ R(t)
```
其中：
- `CBF`: 腦血流量
- `AIF(t)`: 動脈輸入函數
- `R(t)`: 殘留函數 (Residue Function)
- `⊗`: 卷積運算

去卷積的目標是從已知的 C(t) 和 AIF(t) 求解 R(t)，進而計算血流動力學參數。

C(t) 是我們從 CTP 中讀進來的值。AIF(t) 是前面用 huristic 方法自動找出的區域的輸入值。它們都是從一系統 CTP 中算出來的參數。

而 R(t) 就是我們現在需要算出來的結果。它其時代表了生理的模型，也就是腦組織決定了血流要如何經血管到每一個位置的動力學模型。

通過去卷積計算得到重要的血流動力學參數：MTT (平均通過時間)、CBV (腦血容量)、CBF (腦血流量) 和 tMax (最大殘留函數時間)。


### 函數簽名
```python
def runDeconv(imCTC, tIdx, brainMask, aifProps, cSVDThres=0.1, method='bcSVD1', outsideValue=-1):
```

### 功能描述

### 參數說明
- **`imCTC`**: 對比劑時間曲線影像
- **`tIdx`**: 時間索引陣列
- **`brainMask`**: 腦部遮罩 (可以是 SimpleITK Image 或 numpy 陣列)
- **`aifProps`**: AIF 參數 (gamma variate 參數或實際 AIF 曲線)
- **`cSVDThres`**: 循環 SVD 閾值，預設為0.1
- **`method`**: 去卷積方法，預設為 'bcSVD1'，可選 'oSVD', 'bcSVD1', 'bcSVD2'
- **`outsideValue`**: 腦部外區域的填充值，預設為-1



### 主要處理流程

#### 1. 初始化和參數設定
```python
if type(brainMask) is sitk.Image:
    brainMaskNp = sitk.GetArrayFromImage(brainMask)
else:
    brainMaskNp = brainMask

rho, H = 1.05, 0.85  # 腦組織密度和血細胞比容
```
- 處理腦部遮罩的不同輸入格式
- 設定生理常數：腦組織密度 (ρ = 1.05 g/ml) 和血細胞比容 (H = 0.85)

#### 2. AIF 處理
```python
if len(aifProps) != len(tIdx):
    estimAifVal = gv(tIdx, *aifProps)  # 使用 gamma variate 參數生成 AIF
    aifVal = estimAifVal
else:
    aifVal = aifProps  # 直接使用提供的 AIF 曲線
```
- 支援兩種 AIF 輸入格式：gamma variate 參數或實際曲線值

#### 3. 時間間隔計算
```python
deltaT = np.mean(np.diff(tIdx))
```
- 計算平均時間間隔，用於數值積分

#### 4. 卷積矩陣構建
根據選擇的方法構建不同的卷積矩陣：

##### 方法 1: 原始 SVD (oSVD)
```python
if method == 'oSVD':
    G = toeplitz(aifVal, np.zeros(aifVal.shape[0]))
    imCTC_pad = imCTC
```
- 使用標準的 Toeplitz 矩陣
- 不進行零填充

##### 方法 2: 塊循環 SVD 方法 1 (bcSVD1)
```python
elif method == 'bcSVD1':
    # 構建塊循環矩陣的第一列
    colG = np.zeros(2 * tIdx.shape[0])
    colG[0] = aifVal[0]
    colG[tIdx.shape[0]-1] = (aifVal[tIdx.shape[0]-2] + 4 * aifVal[tIdx.shape[0]-1]) / 6
    colG[tIdx.shape[0]] = aifVal[tIdx.shape[0]-1] / 6
    
    for k in range(1, tIdx.shape[0]-1):
        colG[k] = (aifVal[k-1] + 4 * aifVal[k] + aifVal[k+1]) / 6
    
    # 構建第一行
    rowG = np.zeros(2 * tIdx.shape[0])
    rowG[0] = colG[0]
    for k in range(1, 2 * tIdx.shape[0]):
        rowG[k] = colG[2*tIdx.shape[0] - k]
    
    G = toeplitz(colG, rowG)
    imCTC_pad = np.pad(imCTC, [(0, imCTC.shape[0]),] + [(0, 0)]*3)
```
- 使用 Simpson 積分規則構建卷積矩陣
- 對影像進行零填充以減少邊界效應

##### 方法 3: 塊循環 SVD 方法 2 (bcSVD2)
```python
elif method == 'bcSVD2':
    cmat = np.zeros([tIdx.shape[0], tIdx.shape[0]])
    B = np.zeros([tIdx.shape[0], tIdx.shape[0]])
    
    for i in range(tIdx.shape[0]):
        for j in range(tIdx.shape[0]):
            if i == j:
                cmat[i, j] = aifVal[0]
            elif i > j:
                cmat[i, j] = aifVal[(i-j)]
            else:
                B[i, j] = aifVal[tIdx.shape[0] - (j-i)]
    
    G = np.vstack([np.hstack([cmat, B]), np.hstack([B, cmat])])
    imCTC_pad = np.pad(imCTC, [(0, imCTC.shape[0]),] + [(0, 0)]*3)
```
- 另一種塊循環矩陣構建方法
- 同樣進行零填充處理

#### 5. SVD 分解和正則化
```python
U, S, V = np.linalg.svd(G * deltaT)
thres = cSVDThres * np.max(S)
filteredS = 1 / (S + 1e-5)
filteredS[S < thres] = 0
Ginv = V.T @ np.diag(filteredS) @ U.T
```
- 對卷積矩陣進行 SVD 分解
- 使用閾值進行正則化，去除小的奇異值
- 構建偽逆矩陣

#### 6. 殘留函數計算
```python
k = np.abs(np.einsum('ab, bcde->acde', Ginv, imCTC_pad))
k = k[:tIdx.shape[0]]  # 取前半部分（去除填充）
```
- 計算殘留函數 k(t)
- 使用 Einstein 求和約定進行高效矩陣運算

#### 7. 血流動力學參數計算
```python
CBF = H / rho * k.max(axis=0) * 60 * 100        # 腦血流量
CBV = H / rho * k.sum(axis=0) * 100             # 腦血容量
tMax = tIdx[k.argmax(axis=0)]                   # 最大殘留函數時間
MTT = CBV / CBF * 60                            # 平均通過時間
```

各參數的計算公式：
- **CBF**: `(H/ρ) × max(k) × 60 × 100` [ml/100g/min]
- **CBV**: `(H/ρ) × sum(k) × 100` [ml/100g]
- **tMax**: 殘留函數達到最大值的時間 [秒]
- **MTT**: `CBV/CBF × 60` [秒]

#### 8. 遮罩處理
```python
tMax[brainMaskNp == 0] = outsideValue
CBV[brainMaskNp == 0] = outsideValue
CBF[brainMaskNp == 0] = outsideValue
MTT[brainMaskNp == 0] = outsideValue
```
- 將腦部外的區域設為指定填充值

### 返回值
- **`MTT`**: 平均通過時間圖 [秒]
- **`CBV`**: 腦血容量圖 [ml/100g]
- **`CBF`**: 腦血流量圖 [ml/100g/min]
- **`tMax`**: 最大殘留函數時間圖 [秒]

### 方法比較

| 方法 | 優點 | 缺點 | 適用情況 |
|------|------|------|----------|
| oSVD | 計算簡單，速度快 | 邊界效應明顯 | 研究用途，對精度要求不高 |
| bcSVD1 | 減少邊界效應，使用 Simpson 積分 | 計算複雜度中等 | 臨床應用推薦 |
| bcSVD2 | 另一種邊界處理方法 | 計算最複雜 | 特殊情況下使用 |

### 技術特點

#### 優點
1. **多種方法**: 提供三種不同的去卷積方法
2. **正則化處理**: SVD 閾值化避免數值不穩定
3. **邊界效應處理**: bcSVD 方法有效減少邊界偽影
4. **生理參數**: 考慮實際的生理常數

#### 注意事項
1. **SVD 閾值**: 需要根據數據品質調整 cSVDThres
2. **計算複雜度**: bcSVD 方法計算量較大
3. **記憶體需求**: 大影像可能需要大量記憶體
4. **數值穩定性**: 需要適當的正則化參數

### 臨床意義

#### 血流動力學參數的臨床價值
- **CBF**: 反映組織灌注狀態，低值提示缺血
- **CBV**: 反映血管容量，可評估側支循環
- **MTT**: 反映血流通過時間，延長提示灌注障礙
- **tMax**: 類似 TTP，但更穩定，常用於缺血評估

#### 在中風診斷中的應用
- **核心梗死**: CBF 和 CBV 均顯著降低
- **缺血半暗帶**: CBF 降低但 CBV 相對保持
- **側支循環**: CBV 增加，MTT 延長

---

## 總結

`getCTC`、`runTTP`、`runAIF` 和 `runDeconv` 函數構成了完整的 CTP 灌注分析流程：

1. **getCTC**: 提取對比劑濃度變化並進行基線校正
2. **runTTP**: 計算血流動力學參數中的到達峰值時間
3. **runAIF**: 自動識別動脈輸入函數，為後續去卷積分析提供基礎
4. **runDeconv**: 執行去卷積分析，計算核心的血流動力學參數

這個完整的工作流程實現了從原始 DICOM 影像到定量血流動力學參數的全自動分析，在腦中風等疾病的診斷和治療評估中發揮關鍵作用。每個函數都體現了醫學影像處理中的最佳實踐，並提供了足夠的靈活性以適應不同的臨床需求。
