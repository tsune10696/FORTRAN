# ライブラリのインポート
import numpy as np
import matplotlib.pyplot as plt

# データの読み込み
data = np.genfromtxt("output/outputKobe.txt")

# 新しいウィンドウにグラフを作成
plt.figure(1)

# グラフにデータをプロット
plt.plot(data[:,0] * 0.02 , data[:,1], label="responseAcc(m/s/s)", color="red", linestyle="solid", linewidth=1.0)
plt.plot(data[:,0] * 0.02 , data[:,2], label="inputWave(m/s/s)", color="blue", linestyle="solid", linewidth=0.5)

# 最大値／最小値を表示
num = 17996
valuelistAcc = [] #数値を格納する配列
valuelistInput = [] #数値を格納する配列

for i in range(num):
    valueAcc   = data[i,1]
    valueInput = data[i,2]
    valuelistAcc.append(valueAcc)
    valuelistInput.append(valueInput)

print("maxValueAcc=")
print(max(valuelistAcc))
print("minValueAcc=")
print(min(valuelistAcc))
print("maxValueInput=")
print(max(valuelistInput))
print("minValueInput=")
print(min(valuelistInput))

# 軸の調整
plt.xlim([0, 200])  # X軸の範囲を指定
plt.ylim([-20, 20])  # Y軸の範囲を指定

# グラフに情報を表示
plt.title("JMA_KobeWave_NS")
plt.xlabel("time(s)")  # x軸のラベル
plt.ylabel("acc(m/s/s)")  # y軸のラベル
plt.legend(loc="lower right")  # 汎例
plt.grid(True)  # グリッド線

# グラフの描画を実行
plt.savefig("output/equationOfMotion.png")
plt.show()
