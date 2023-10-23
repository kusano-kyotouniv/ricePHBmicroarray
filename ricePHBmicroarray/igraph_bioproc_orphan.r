library(igraph)		# R ライブラリ "igraph" をロードする
data1 <- read.table(file="bigraph_bioproc_orphan.txt", sep = "\t",header=T, quote = '"')		# bigraphフォーマット？のデータをロードする
data2 <- read.table(file="vertex_bioproc_orphan.txt", sep = "\t", header=T, quote = '"')

#subdat <- subset(data, freq>0)		# ヘッダーが"freq" の列が 0.4 より大きいものだけ選抜
net<-graph.data.frame(data1, vertices=data2, directed = T)		# plot関数の引数に変換する。directed=FALSE で矢印なし

pdf("image_bioproc_orphan.pdf",width=40,height=40)		# 出力を画像ファイル(.pdf)にする。画像サイズは文字のサイズに影響する
plot(net,vertex.size = 3,#E(net)$freq*2,	# ○のサイズ。データ依存にできるはず E(net)$freq*4
		 vertex.color = V(net)$color,	# "white", 		# ○の塗り色
		 vertex.label = V(net)$name,
		 vertex.label.font = 2, 		# ラベル文字をboldモード(2番)にする。3だとイタリックとか
		 vertex.label.cex = 1.5, 		# ラベル文字のサイズ
		 vertex.label.color = "blue", 	# ラベル文字の色 "blue"
		 vertex.label.dist = 0, 		# ラベル文字の表示位置を○からちょっとずらすモード（１番）
		 edge.arrow.size = 0.8, 		# 矢印の矢の部分のサイズ調整
		 edge.width=E(net)$freq, 				# 線の太さをデータ依存にする。"freq" は元データのヘッダーの記号
		 layout=layout.fruchterman.reingold )	# レイアウトアルゴリズムを選ぶ。いろいろあるらしい
dev.off()								# pdf関数を終了。たぶんこれないとファイルが出ない
