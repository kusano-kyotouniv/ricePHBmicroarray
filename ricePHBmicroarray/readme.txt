Shimada et al., イネカルスでpolyhydroxybutanoate を作ろうとした実験の
マイクロアレイ解析部分の情報解析に関するMethod。
各Figureに示される実験結果は
以下1-5のメソッドをこの順番で実行することで得られた。

Fig4A, 階層クラスタリング
	1 -> 2 

Fig6 と Fig4B-D, 散布図
	1 -> 3　-> マウスで操作

Fig5, 共起ネットワーク
	1 -> 4 -> 5

なお、共起ネットワーク図の線分の長さと角度は
計算するごとに異なるが、接続される頂点は同じである。
これは使用したレイアウトがランダム関数を含むためである。


＊＊＊
以下、1-5に各メソッドの説明を記す。


1. make_GOenrich_gene.pl

アノテーション情報ファイルとマイクロアレイデータをロードして
GOエンリッチメント解析の結果を出力するperlスクリプト。

ロードするファイル
	IRGSP1.txt	アノテーション情報ファイル。RAPの22年9月版
	raw1.txt	マイクロアレイの実験結果ファイル
	raw2.txt	マイクロアレイの実験結果ファイル
	raw3.txt	マイクロアレイの実験結果ファイル
	raw4.txt	マイクロアレイの実験結果ファイル
	※注意：実験結果ファイルには今回使っていない実験のデータも含まれる。

生成するファイル
	tamakiarray_hclust.txt
		階層クラスタリング用のデータ。
		hclust.r でロードされる。

	tamakiarray_gene_RNAseqViewer2.txt
		scatterplot用の遺伝子発現データ。
		volcano plot のP値を含む。
		RNAseqViewer6_for_mac.pl でRSEMdataとしてロードする

	GOenrich_molfunc.txt
	GOenrich_celcomp.txt
	GOenrich_bioproc.txt
		network図を生成するためのデータ。
		Molecular Function, 
		Biological Process,
		Cellular Component の３カテゴリ。



2. hclust.r

階層クラスタリングを計算するRスクリプト。

ロードするファイル
	tamakiarray_hclust.txt	1.で生成したファイル

生成するファイル
	hclust.pdf		階層クラスタリング結果の画像


3. RNAseqViewer6_for_mac.pl

トランスクリプトームデータを見るためのGUIソフトウェア。
遺伝子発現量のデータを散布図で表示する機能がある。
散布図は表示されると同時にPNG形式の画像ファイルとして出力される。

ロードするファイル
	RNAseqViewer_settings.txt			散布図に使うデータがどれかを指定する内部ファイル	
	tamakiarray_gene_RNAseqViewer2.txt	1.で生成したファイル。
		[Load RSEM] ボタンをマウスでクリックしてこのファイルを指定すると
		散布図の画像ファイルが生成する。

生成するファイル
	RNAseqViewer_scatplot.png			散布図の画像ファイル


4. make_bigraph2.pl

ネットワークグラフ図を生成するための共起情報（bigraph）ファイルを生成する

ロードするファイル
	transcripts.gff			アノテーション情報ファイル。RAPの22年9月版
	GOenrich_molfunc.txt	1.で生成したファイル。
	GOenrich_celcomp.txt	1.で生成したファイル。
	GOenrich_bioproc.txt	1.で生成したファイル。

生成するファイル
	bigraph_molfunc.txt				メインの共起情報データ、Molecular Function 枠		
	bigraph_molfunc_rare.txt		１個しか共起GOの無いもの、Molecular Function 枠
	bigraph_molfunc_orphan.txt		共起GOが１つも無いもの、Molecular Function 枠
	bigraph_bioproc.txt				Biological Process 枠の同セット。以下略
	bigraph_bioproc_rare.txt		
	bigraph_bioproc_orphan.txt
	bigraph_celcomp.txt				Cellular Component 枠の同セット。以下略
	bigraph_celcomp_rare.txt
	bigraph_celcomp_orphan.txt

	vertex_molfunc.txt				各bigraphデータの頂点色指定データ。以下略
	vertex_molfunc_rare.txt
	vertex_molfunc_orphan.txt
	vertex_bioproc.txt
	vertex_bioproc_rare.txt
	vertex_bioproc_orphan.txt
	vertex_celcomp.txt
	vertex_celcomp_rare.txt
	vertex_celcomp_orphan.txt


5. run_igraph.sh

bigraphデータをロードしてネットワークグラフ図を生成する
Rスクリプトを実行するためのシェルスクリプト

ロードするファイル
	igraph_bioproc_orphan.r		4.で生成したファイルをロードしてpdf画像を出力する
	igraph_bioproc_rare.r		以下略
	igraph_bioproc.r
	igraph_celcomp_orphan.r
	igraph_celcomp_rare.r
	igraph_celcomp.r
	igraph_molfunc_orphan.r
	igraph_molfunc_rare.r
	igraph_molfunc.r
	4. で生成したファイル	

生成するファイル
	image_bioproc_orphan.pdf	ネットワークグラフ図のpdfファイル。
	image_bioproc_rare.pdf
	image_bioproc.pdf
	image_celcomp_orphan.pdf
	image_celcomp_rare.pdf
	image_celcomp.pdf
	image_molfunc_orphan.pdf
	image_molfunc_rare.pdf
	image_molfunc.pdf

	
