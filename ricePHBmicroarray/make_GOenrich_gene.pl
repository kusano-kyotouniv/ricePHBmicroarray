# use Math::Trig 'pi';		# P値の計算に円周率が出てきた。でも要らなかった。
$| = 1;

# ４個めのスクリプト。
# アノテーションを新しいものにしてみる。割と使えるスポットが増えた。GO termもちょっと増えた？
# RNAseqViewerで表示する用のデータも必要になったのでついでに出力する。
# 遺伝子単位のデータもついでに出力する。めんどくさい
# でも、これで古い方のデータもRNAseqViewerで表示できるようになるな。便利

# アノテーション情報ファイルをロードする。
open ANN, 'IRGSP1.txt';		# アノテーション情報のファイル。ロードしやすくトリミングずみ
#open ANN, 'RAP3_trunc.txt';		# アノテーション情報のファイル。ロードしやすくトリミングずみ
my %an;
my $ans=0;
while(my $line = <ANN>){
	$line =~ s/\r|\n|\r\n//g;			# 改行コードが入ってるみたい。windowsのコードか心配なので取る
	my @s=split(/\t/,$line);
	$an{$s[0]}=$s[1];			# ゼロ番目がIDで、１番めがデータ本体らしい。２列なのかな。
	$ans++;
}
close ANN;
print 'RAP3 annotation entries: '.$ans."\n";

# マイクロアレイデータをロードする。まずヘッダーをスキップ
print 'loading microarray data...';
open IN1, 'raw1.txt';			# マイクロアレイ結果の元データファイル
open IN2, 'raw2.txt';
open IN3, 'raw3.txt';
open IN4, 'raw4.txt';
for(my $i=0;$i<21;$i++){		# 最初の21行はヘッダーらしい。
	my $line1=<IN1>;
	my $line2=<IN2>;
	my $line3=<IN3>;
	my $line4=<IN4>;
}

# マイクロアレイデータをロードする。本体。
my $blankspot13 = 'DarkCorner|art|GE_BrightCorner|random|genome|NegativeControl|EQC';
my $blankspot4 = 'POsControl';
my $rap_hits=0;
my $count=0;
my $lines=0;
my @gid;		# $gid[genes]	遺伝子ID。Os##g######
my @acce;		# $acce[genes]	アクセッション記号。遺伝子IDに追加でなんか描いてある
my @flcd;		# $flcd[genes]	FLcDNA のIDとかかいてあるやつ
my @anno;		# $anno[genes]	簡易アノテーション。RAPの方は後で足す。
my @fold1;		# 発現変動。コントロールはraw4 のgreen
my @fold2;		# $fold_[genes][spots] まず割り算の結果を配列で撮っておく。
my @fold3;		# あとで割り算の結果を足し算で積算する。最後に平均にする。
my @fold4;
my @spots;		# $spots[genes]	遺伝子ごとのスポット数。平均とるときの分母。
my %index;		# $index{id} = genes
my $genes=0;
my @sign1;		# $sign[genes][spots]
my @sign2;		# 1=raw1の赤、2=raw2の赤、3=raw3の赤、4=raw4の赤、5=raw4の緑。
my @sign3;
my @sign4;
my @sign5;		# 8=raw4の緑。raw1の緑も使ってみたいので@sign5-8とする。とりあえずgeneの分だけ
my @sign6;
my @sign7;
my @sign8;
open OUT, '>tamakiarray_spot_RNAseqViewer.txt';
print OUT 'Contig ID'."\t";
print OUT 'Contig length'."\t";
print OUT 'SpkBC1 (raw1) TPM'."\t";		# TPM はRNAseqViewerでデータ列を認識するための目印
print OUT 'SpkBC2 (raw2) TPM'."\t";
print OUT '35S-C1 (raw3) TPM'."\t";
print OUT '35S-C2 (raw4) TPM'."\t";
print OUT 'empty vector control (raw1) TPM'."\t";
print OUT 'empty vector control (raw2) TPM'."\t";
print OUT 'empty vector control (raw3) TPM'."\t";
print OUT 'empty vector control (raw4) TPM'."\t";
print OUT 'PhaC expressing (AU) TPM'."\t";
print OUT 'empty vector control (AU) TPM'."\t";
print OUT 'annotation'."\t";
print OUT 'IRGSP1 annotation'."\n";

open OUT2, '>tamakiarray_hclust.txt';
print OUT2 'spot no'."\t";
print OUT2 'proSPK-BC1'."\t";
print OUT2 'proSPK-BC2'."\t";
print OUT2 'pro35S-C1'."\t";
print OUT2 'pro35S-C2'."\t";
print OUT2 'empty vector1'."\t";
print OUT2 'empty vector2'."\t";
print OUT2 'empty vector3'."\t";
print OUT2 'empty vector4'."\n";

WHL:while(	my $line1 = <IN1>){			# データを一行ずつロードしていく
			my $line2 = <IN2>;
			my $line3 = <IN3>;
			my $line4 = <IN4>;
		
	my @s1 = split(/\t/,$line1);		# データはタブ区切りらしい。助かる。
	my @s2 = split(/\t/,$line2);
	my @s3 = split(/\t/,$line3);
	my @s4 = split(/\t/,$line4);
	
	my $spotnum = $s1[1];		# １番はスポット番号
	my $accession = $s1[4];		# ４番がアクセッション番号
	my $flcdna = $s1[13];		# １３番がFLcDNAのアクセッション番号
	my $annotation = $s1[14];	# １４番はアノテーション情報。
	
	$lines++;
	if($accession =~ /$blankspot4/){next;}		# 使いようのないスポットは処理しない
	if($accession eq ''){next;}
	if($flcdna =~ /$blankspot13/){next;}
	
	my $signal_g1 = $s1[24]/1000;	# ２４と２５番は緑と赤	RNAseqViewer仕様にするため/1000する。ミリ輝度値
	my $signal_r1 = $s1[25]/1000;
	my $signal_g2 = $s2[24]/1000;	# ２４と２５番は緑と赤
	my $signal_r2 = $s2[25]/1000;
	my $signal_g3 = $s3[24]/1000;	# ２４と２５番は緑と赤
	my $signal_r3 = $s3[25]/1000;
	my $signal_g4 = $s4[24]/1000;	# ２４と２５番は緑と赤
	my $signal_r4 = $s4[25]/1000;
	
	# 階層クラスタリング用のデータ。最低値ゼロでもOKなので最低値を設定する前に記録する
	print OUT2 $spotnum."\t";
	print OUT2 $signal_r1."\t";
	print OUT2 $signal_r2."\t";
	print OUT2 $signal_r3."\t";
	print OUT2 $signal_r4."\t";
	print OUT2 $signal_g1."\t";
	print OUT2 $signal_g2."\t";
	print OUT2 $signal_g3."\t";
	print OUT2 $signal_g4."\n";

	# RNAseqViewer で Log軸を使いたいので最低値を設定する。
	if($signal_g1<0.01){$signal_g1=0.01}		# 小さすぎるシグナル値は最低値のリミット値を設定して上げる
	if($signal_g2<0.01){$signal_g2=0.01}		# 生データでシグナル値10以下は暴れるので、ミリ輝度値の0.01を最低値とする
	if($signal_g3<0.01){$signal_g3=0.01}
	if($signal_g4<0.01){$signal_g4=0.01}
	if($signal_r1<0.01){$signal_r1=0.01}
	if($signal_r2<0.01){$signal_r2=0.01}
	if($signal_r3<0.01){$signal_r3=0.01}
	if($signal_r4<0.01){$signal_r4=0.01}
	
#	my $id = substr($accession,0,12);
	my @ii = split(/\|/,$accession);
	my $id = $ii[0];
	my $control_average = ($signal_g1+$signal_g2+$signal_g3+$signal_g4)/4;
	my $experiment_average = ($signal_r1+$signal_r2+$signal_r3+$signal_r4)/4;
	if(exists($index{$id})){			# まず遺伝子indexを探す。既存だったら足し算。
		my $g = $index{$id};	# hashは遅そうなので遺伝子のindexはスカラに一旦取ってみる
		$fold1[$g][$spots[$g]] = $signal_r1/$control_average;		# まず二次元配列にデータを取る。平均はあとで計算する。
		$fold2[$g][$spots[$g]] = $signal_r2/$control_average;
		$fold3[$g][$spots[$g]] = $signal_r3/$control_average;
		$fold4[$g][$spots[$g]] = $signal_r4/$control_average;
		$sign1[$g][$spots[$g]] = $signal_r1;
		$sign2[$g][$spots[$g]] = $signal_r2;
		$sign3[$g][$spots[$g]] = $signal_r3;
		$sign4[$g][$spots[$g]] = $signal_r4;
		$sign5[$g][$spots[$g]] = $signal_g1;
		$sign6[$g][$spots[$g]] = $signal_g2;
		$sign7[$g][$spots[$g]] = $signal_g3;
		$sign8[$g][$spots[$g]] = $signal_g4;
		$spots[$g]++;
	}else{
		$index{$id} = $genes;			# 新規だったら新規登録。
		$fold1[$genes][0] = $signal_r1/$control_average;
		$fold2[$genes][0] = $signal_r2/$control_average;
		$fold3[$genes][0] = $signal_r3/$control_average;
		$fold4[$genes][0] = $signal_r4/$control_average;
		$sign1[$genes][0] = $signal_r1;
		$sign2[$genes][0] = $signal_r2;
		$sign3[$genes][0] = $signal_r3;
		$sign4[$genes][0] = $signal_r4;
		$sign5[$genes][0] = $signal_g1;
		$sign6[$genes][0] = $signal_g2;
		$sign7[$genes][0] = $signal_g3;
		$sign8[$genes][0] = $signal_g4;
		$gid[$genes]=$id;	# 今回使うのはこっち。
		$acce[$genes]=$id;	# こっちは一応とっとく。使わない予定。
		$flcd[$genes]=$flcdna;
		$anno[$genes]=$annotation;
		$spots[$genes]=1;
		$genes++;
	}
	
	# RNAseqViewer用のデータのためにアノテーション情報を照合
	$rap_annotation='';
	if(exists($an{$id})){$rap_annotation=$an{$id};}
	
	# RNAseqViewer用のデータを出力。spot単位の方
	print OUT $id."\t";
	print OUT '1'."\t";		# ２列めは色に使う値。本来はcontig length。どうしようかな。ゼロでない数が入用。とりあえず１で。
	print OUT $signal_r1."\t";
	print OUT $signal_r2."\t";
	print OUT $signal_r3."\t";
	print OUT $signal_r4."\t";
	print OUT $signal_g1."\t";
	print OUT $signal_g2."\t";
	print OUT $signal_g3."\t";
	print OUT $signal_g4."\t";
	print OUT $experiment_average."\t";
	print OUT $control_average."\t";
	print OUT $annotation."\t";
	print OUT $rap_annotation."\n";

	$count++;
	if($count % 1000 == 0){print '.';}
}
close OUT2;
close OUT;
close IN4;
close IN3;
close IN2;
close IN1;
print 'done'."\n";

print 'total spots found in the microarray: '.$lines."\n";
print 'experimental spots in the microarray: '.$count."\n";
print 'total genes found including ncRNAs: '.$genes."\n";

### ロード終了。まずfoldchange の平均を計算する
my @ave_fc1;		# $ave_fc_[genes]
my @ave_fc2;
my @ave_fc3;
my @ave_fc4;
for(my $g=0;$g<$genes;$g++){
	for(my $p=0;$p<$spots[$g];$p++){
		$ave_fc1[$g] += $fold1[$g][$p];
		$ave_fc2[$g] += $fold2[$g][$p];
		$ave_fc3[$g] += $fold3[$g][$p];
		$ave_fc4[$g] += $fold4[$g][$p];
	}
	$ave_fc1[$g] = $ave_fc1[$g]/$spots[$g];
	$ave_fc2[$g] = $ave_fc2[$g]/$spots[$g];
	$ave_fc3[$g] = $ave_fc3[$g]/$spots[$g];
	$ave_fc4[$g] = $ave_fc4[$g]/$spots[$g];
}

### RAP3 と照合する。
my @rapa;	# $rapa[genes]
for(my $g=0;$g<$genes;$g++){
	if(exists($an{$gid[$g]})){$rapa[$g] = $an{$gid[$g]};$rap_hits++;}
}
print 'genes with RAP3 annotation hit: '.$rap_hits."\n";

# foldchange値のデータを出力
open OUT, '>tamakiarray_gene_foldchange_excel.txt';

print OUT 'Accession number'."\t";
print OUT 'spots'."\t";
print OUT 'SpkBC1 foldchange'."\t";
print OUT 'SpkBC2 foldchange'."\t";
print OUT '35SC1 foldchange'."\t";
print OUT '35SC2 foldchange'."\t";
print OUT 'SpkBC average'."\t";
print OUT '35SC average'."\t";
print OUT 'Annotation'."\t";
print OUT 'GO Molecular Function'."\t";
print OUT 'GO Cellular Component'."\t";
print OUT 'GO Biological Process'."\t";
print OUT 'GO Molecular Function'."\t";
print OUT 'GO Cellular Component'."\t";
print OUT 'GO Biological Process'."\t";
print OUT 'IRGSP1 annotation'."\t";
print OUT "\n";

for(my $g=0;$g<$genes;$g++){
	my $got_molfunc = '';		# GOterm を取り出して表示したい
	my $got_celcomp = '';
	my $got_bioproc = '';
	my $goi_molfunc = '';
	my $goi_celcomp = '';
	my $goi_bioproc = '';
	my $go_molfunc = '';
	my $go_celcomp = '';
	my $go_bioproc = '';
	
	if($rapa[$g] =~ /GO\=.+?\;/){
		my $text = $&;
		while($text =~ m/Molecular Function\:.+?GO\:[0-9]{7,7}\)/g){
			my $go = $&;
			$goi_molfunc = $goi_molfunc.'|'.substr($go,-11,-1);	# GO:####### を並べる
			$got_molfunc = $got_molfunc.'|'.substr($go,20,-13);	# GO term を並べる
		}
		while($text =~ m/Cellular Component\:.+?GO\:[0-9]{7,7}\)/g){
			my $go = $&;
			$goi_celcomp = $goi_celcomp.'|'.substr($go,-11,-1);	# GO:####### を並べる
			$got_celcomp = $got_celcomp.'|'.substr($go,20,-13);	# GO term を並べる
		}
		while($text =~ m/Biological Process\:.+?GO\:[0-9]{7,7}\)/g){
			my $go = $&;
			$goi_bioproc = $goi_bioproc.'|'.substr($go,-11,-1);	# GO:####### を並べる
			$got_bioproc = $got_bioproc.'|'.substr($go,20,-13);	# GO term を並べる
		}
		$goi_molfunc = substr($goi_molfunc,1);
		$goi_celcomp = substr($goi_celcomp,1);
		$goi_bioproc = substr($goi_bioproc,1);
		$got_molfunc = substr($got_molfunc,1);
		$got_celcomp = substr($got_celcomp,1);
		$got_bioproc = substr($got_bioproc,1);
		if($got_molfunc ne ''){$go_molfunc = $got_molfunc.'; '.$goi_molfunc}
		if($got_celcomp ne ''){$go_celcomp = $got_celcomp.'; '.$goi_celcomp}
		if($got_bioproc ne ''){$go_bioproc = $got_bioproc.'; '.$goi_bioproc}
	}

	print OUT $gid[$g]."\t";
	print OUT $spots[$g]."\t";
	print OUT $ave_fc1[$g]."\t";
	print OUT $ave_fc2[$g]."\t";
	print OUT $ave_fc3[$g]."\t";
	print OUT $ave_fc4[$g]."\t";
	print OUT (($ave_fc1[$g]+$ave_fc2[$g])/2)."\t";
	print OUT (($ave_fc3[$g]+$ave_fc4[$g])/2)."\t";
	print OUT $anno[$g]."\t";
	print OUT $got_molfunc."\t";
	print OUT $got_celcomp."\t";
	print OUT $got_bioproc."\t";
	print OUT $goi_molfunc."\t";
	print OUT $goi_celcomp."\t";
	print OUT $goi_bioproc."\t";
	print OUT $rapa[$g]."\n";
}
close OUT;

# 輝度値のデータを出力する。散布図表示用。control１個だったときの旧版
open OUT, '>tamakiarray_gene_RNAseqViewer.txt';
print OUT 'Contig ID'."\t";
print OUT 'Contig length'."\t";
print OUT 'SpkBC1 (raw1) TPM'."\t";
print OUT 'SpkBC2 (raw2) TPM'."\t";
print OUT '35S-C1 (raw3) TPM'."\t";
print OUT '35S-C1 (raw4) TPM'."\t";
print OUT 'empty vector control (raw4) TPM'."\t";
print OUT 'annotation'."\t";
print OUT 'IRGSP1 annotation'."\n";

# 論文用。scatplotするためのやつ。volcano plot にしたいけど流石に時間がないか
open OUT2, '>tamakiarray_gene_RNAseqViewer2.txt';
print OUT2 'Gene ID'."\t";
print OUT2 'responsible spots'."\t";
print OUT2 'PhaC expressing (AU) expression'."\t";
print OUT2 'empty vector control (AU) expression'."\t";
print OUT2 'PhaC/control foldchange expression'."\t";
print OUT2 'P-value expression'."\t";
#print OUT2 'empty vector control (AU) expression'."\t";
print OUT2 'SpkBC expressing (AU) expression'."\t";
print OUT2 '35SC expressing (AU) expression'."\t";
print OUT2 'SpkBC expressing1 (AU) expression'."\t";
print OUT2 'SpkBC expressing2 (AU) expression'."\t";
print OUT2 '35SC expressing1 (AU) expression'."\t";
print OUT2 '35SC expressing2 (AU) expression'."\t";
print OUT2 'empty vector control1 (AU) expression'."\t";
print OUT2 'empty vector control2 (AU) expression'."\t";
print OUT2 'empty vector control3 (AU) expression'."\t";
print OUT2 'empty vector control4 (AU) expression'."\t";
print OUT2 'original annotation'."\t";
print OUT2 'IRGSP1 latest annotation'."\n";

# my $err=0;	#分散の計算２種類が本当に同じか調べるカウンター。合ってたのでもう要らない
my @t_value;
my @p_value;
my $s_max=0;		# 全部の輝度値が等しいと分散がゼロになって割り算できなくなるので最小値を調べる
my $s_min=1000;		# ついでに最大値も調べる

my $t_dist_factor = 2.5*1.5*0.5/2/sqrt(6);	# 自由度６のときのt分布の式の係数。サンプル数-2（２群の場合の自由度）から一義的に定まる
my $t_step = 0.1;		# t分布の計算で積分をするときの精度。台形の面積を足していく。台形の高さにあたる値
my $t_end = 30;			# t分布の端を雑に決める。本当は無限大にしたいが無理なので。データ中でt値の最大が14くらいのようなので、その倍くらい。
my $t_area_total=0;		# P値を計算するときの分母。P値の最大値は100%なので。t分布をゼロから端まで積分。定数なので先に計算しておく。
for(my $i=0;$i<$t_end;$i+=$t_step){
	my $left = $t_dist_factor*(1+($i**2)/6)**(-3.5);
	my $right = $t_dist_factor*(1+(($i+$t_step)**2)/6)**(-3.5);
	my $area = ($left+$right)*$t_step/2;
	$t_area_total += $area;
}
print 't_area_total: '.$t_area_total."\n";

for(my $g=0;$g<$genes;$g++){
	print OUT $gid[$g]."\t";
	print OUT $spots[$g]."\t";
	my $intensity1=0;		# 輝度値は遺伝子に所属するスポットの平均値
	my $intensity2=0;
	my $intensity3=0;
	my $intensity4=0;
	my $intensity5=0;		# 5は一枚scatplot用のデータ出力の計算のためのもの
	my $intensity6=0;		# 5は一枚scatplot用のデータ出力の計算のためのもの
	my $intensity7=0;		# 5は一枚scatplot用のデータ出力の計算のためのもの
	my $intensity8=0;
	for(my $p=0;$p<$spots[$g];$p++){
		$intensity1 += $sign1[$g][$p];
		$intensity2 += $sign2[$g][$p];
		$intensity3 += $sign3[$g][$p];
		$intensity4 += $sign4[$g][$p];
		$intensity5 += $sign5[$g][$p];
		$intensity6 += $sign6[$g][$p];
		$intensity7 += $sign7[$g][$p];
		$intensity8 += $sign8[$g][$p];
	}
	$intensity1 = $intensity1/$spots[$g];
	$intensity2 = $intensity2/$spots[$g];
	$intensity3 = $intensity3/$spots[$g];
	$intensity4 = $intensity4/$spots[$g];
	$intensity5 = $intensity5/$spots[$g];
	$intensity6 = $intensity6/$spots[$g];
	$intensity7 = $intensity7/$spots[$g];
	$intensity8 = $intensity8/$spots[$g];
	
	print OUT $intensity1."\t";
	print OUT $intensity2."\t";
	print OUT $intensity3."\t";
	print OUT $intensity4."\t";
	print OUT $intensity8."\t";
	print OUT $anno[$g]."\t";
	print OUT $rapa[$g]."\n";

	print OUT2 $gid[$g]."\t";
	print OUT2 $spots[$g]."\t";
	my $intensity_expr=0;		# 輝度値は遺伝子に所属するスポットの平均値
	my $intensity_expr_SPKBC=0;
	my $intensity_expr_35S_C=0;
	my $intensity_cont=0;
	my $intensity_cont_58=0;
	for(my $p=0;$p<$spots[$g];$p++){
		$intensity_expr += $sign1[$g][$p];
		$intensity_expr += $sign2[$g][$p];
		$intensity_expr += $sign3[$g][$p];
		$intensity_expr += $sign4[$g][$p];
		$intensity_expr_SPKBC += $sign1[$g][$p];
		$intensity_expr_SPKBC += $sign2[$g][$p];
		$intensity_expr_35S_C += $sign3[$g][$p];
		$intensity_expr_35S_C += $sign4[$g][$p];
		$intensity_cont += $sign5[$g][$p];
		$intensity_cont += $sign6[$g][$p];
		$intensity_cont += $sign7[$g][$p];
		$intensity_cont += $sign8[$g][$p];
		$intensity_cont_58 += $sign5[$g][$p];
		$intensity_cont_58 += $sign8[$g][$p];
	}
	$intensity_expr = $intensity_expr/$spots[$g]/4;
	$intensity_expr_SPKBC = $intensity_expr_SPKBC/$spots[$g]/2;
	$intensity_expr_35S_C = $intensity_expr_35S_C/$spots[$g]/2;
	$intensity_cont = $intensity_cont/$spots[$g]/4;
	$intensity_cont_58 = $intensity_cont_58/$spots[$g]/2;
	
	$foldchange = $intensity_expr/$intensity_cont;

	# ここで分散とP値を計算してみたい。２群間のt検定で対応のある場合は、分散は実験区とコントロールのデータを一緒くたして算出する。
	my $sw=1;
	my $p_val=0;
	if($sw==1){
		my $total_average = ($intensity_expr + $intensity_cont )/2;		# 平均値を使うので計算しておく

		my $s = 0; # 分散 = sigma(sample_value^2)/samples - average^2	# こっちのほうが早いらしい。
		for(my $p=0;$p<$spots[$g];$p++){
			$s += $sign1[$g][$p]**2;
			$s += $sign2[$g][$p]**2;
			$s += $sign3[$g][$p]**2;
			$s += $sign4[$g][$p]**2;
			$s += $sign5[$g][$p]**2;
			$s += $sign6[$g][$p]**2;
			$s += $sign7[$g][$p]**2;
			$s += $sign8[$g][$p]**2;
		}
		$s = $s/$spots[$g]/8 - $total_average**2;
		if($s_max < $s){$s_max = $s;}
		if($s_min > $s){$s_min = $s;}	# $s の最低値を調べたい。ゼロのはず・・

		# ゼロで割ったエラー防止。最低値を割り当てたい。$s = sum{(0.01-0.01)^2}/8 = 0 の場合のはず。なので、0.01 の２乗の 0.0001 を割り当てる。
		if($s<0.0001){$s = 0.0001;}	

		# 次にt値を計算する。
		my $t = abs($intensity_expr - $intensity_cont)/sqrt($s/$spots[$g]/4);	# ゼロ割エラーがここで出る。
		$t_value[$g] = $t;

		# そしてP値に変換する。
		$p_val = p_value($t);
		$p_value[$g] = $p_val;

		# 計算合ってるか見てみる
		if($g % 2000 == 0){				# 計算合ってるか調べる。
			my $ii = int($intensity_expr*1000)/1000;
			my $ss = int($s*10000)/10000;
			my $tt = int($t*100)/100;
			my $pp = int($p_val*10000)/10000;
			my $ppp = int(-log($pp)/log(10)*100)/100;
			print $g."\t".$ii."\t".$ss."\t".$tt."\t".$pp."\t".$ppp."\n";
		}
	}


	print OUT2 $intensity_expr."\t";
	print OUT2 $intensity_cont."\t";
	print OUT2 $foldchange."\t";
	print OUT2 $p_val."\t";
	print OUT2 $intensity_expr_SPKBC."\t";
	print OUT2 $intensity_expr_35S_C."\t";
#	print OUT2 $intensity_cont_58."\t";
	print OUT2 $intensity1."\t";
	print OUT2 $intensity2."\t";
	print OUT2 $intensity3."\t";
	print OUT2 $intensity4."\t";
	print OUT2 $intensity5."\t";
	print OUT2 $intensity6."\t";
	print OUT2 $intensity7."\t";
	print OUT2 $intensity8."\t";
	print OUT2 $anno[$g]."\t";
	print OUT2 $rapa[$g]."\n";


}
close OUT2;
close OUT;

# P値の計算ができたかチェク
print 's_max: '.$s_max."\n";
print 's_min: '.$s_min."\n";

my @t_hist;
my $t_max=0;
for(my $g=0;$g<$genes;$g++){
	my $tt = int($t_value[$g]);
	$t_hist[$tt]++;
	if($t_max<$tt){$t_max=$tt;}
}
print 't-value distribution: '."\n";
for(my $i=0;$i<$t_max+1;$i++){print $i."\t".$t_hist[$i]."\n";}
print "\n";

my @p_hist;
my $p_max=0;
for(my $g=0;$g<$genes;$g++){
	my $pp = int((-log($p_value[$g]))/log(10));
	$p_hist[$pp]++;
	if($p_max<$pp){$p_max=$pp;}
}
print 'p-value distribution (-log10): '."\n";
for(my $i=0;$i<$p_max+1;$i++){print $i."\t".$p_hist[$i]."\n";}


##### GOを取り出す。とりあえずfoldchenge を掛け算で積算。
my %goterm_molfunc;
my %goterm_celcomp;
my %goterm_bioproc;
my @goid_molfunc;
my @goid_celcomp;
my @goid_bioproc;
my %goindex_molfunc;
my %goindex_celcomp;
my %goindex_bioproc;
my @gogenes_molfunc;
my @gogenes_celcomp;
my @gogenes_bioproc;
my $gos_molfunc=0;
my $gos_celcomp=0;
my $gos_bioproc=0;
my @goenrich1_molfunc;
my @goenrich1_celcomp;
my @goenrich1_bioproc;
my @goenrich2_molfunc;
my @goenrich2_celcomp;
my @goenrich2_bioproc;
my @goenrich3_molfunc;
my @goenrich3_celcomp;
my @goenrich3_bioproc;
my @goenrich4_molfunc;
my @goenrich4_celcomp;
my @goenrich4_bioproc;
my $gene_with_gos=0;
for(my $g=0;$g<$genes;$g++){
	if($rapa[$g] =~ /GO\=.+?\;/){
		$gene_with_gos++;
		my $text = $&;
		while($text =~ m/Molecular Function\:.+?GO\:[0-9]{7,7}\)/g){
			my $go = $&;
			my $goid = substr($go,-11,-1);		# GO:####### を取り出す。
			if(exists($goterm_molfunc{$goid})){
				my $index = $goindex_molfunc{$goid};	# hash は遅そうなので取り出すのは１回にす	る
				$goenrich1_molfunc[$index] *= $ave_fc1[$g];
				$goenrich2_molfunc[$index] *= $ave_fc2[$g];			# foldchange 値は掛け算で積算する
				$goenrich3_molfunc[$index] *= $ave_fc3[$g];
				$goenrich4_molfunc[$index] *= $ave_fc4[$g];
				$gogenes_molfunc[$index]++;
			}else{
				$goindex_molfunc{$goid} = $gos_molfunc;			# 配列のindexを登録
				$goid_molfunc[$gos_molfunc] = $goid;			# GO番号 を登録
				$goterm_molfunc{$goid} = substr($go,20,-13);	# GO term を取り出す。
				$gogenes_molfunc[$gos_molfunc]=1;
				$goenrich1_molfunc[$gos_molfunc] = $ave_fc1[$g];
				$goenrich2_molfunc[$gos_molfunc] = $ave_fc2[$g];	# foldchange 値を４実験区分
				$goenrich3_molfunc[$gos_molfunc] = $ave_fc3[$g];
				$goenrich4_molfunc[$gos_molfunc] = $ave_fc4[$g];
				$gos_molfunc++;
			}
		}
		while($text =~ m/Cellular Component\:.+?GO\:[0-9]{7,7}\)/g){
			my $go = $&;
			my $goid = substr($go,-11,-1);		# GO:####### を取り出す。
			if(exists($goterm_celcomp{$goid})){
				my $index = $goindex_celcomp{$goid};	# hash は遅そうなので取り出すのは１回にする
				$goenrich1_celcomp[$index] *= $ave_fc1[$g];
				$goenrich2_celcomp[$index] *= $ave_fc2[$g];			# foldchange 値は掛け算で積算する
				$goenrich3_celcomp[$index] *= $ave_fc3[$g];
				$goenrich4_celcomp[$index] *= $ave_fc4[$g];
				$gogenes_celcomp[$index]++;
			}else{
				$goindex_celcomp{$goid} = $gos_celcomp;			# 配列のindexを登録
				$goid_celcomp[$gos_celcomp] = $goid;			# GO番号 を登録
				$goterm_celcomp{$goid} = substr($go,20,-13);	# GO term を取り出す。
				$gogenes_celcomp[$gos_celcomp]=1;
				$goenrich1_celcomp[$gos_celcomp] = $ave_fc1[$g];
				$goenrich2_celcomp[$gos_celcomp] = $ave_fc2[$g];	# foldchange 値を４実験区分
				$goenrich3_celcomp[$gos_celcomp] = $ave_fc3[$g];
				$goenrich4_celcomp[$gos_celcomp] = $ave_fc4[$g];
				$gos_celcomp++;
			}
		}
		while($text =~ m/Biological Process\:.+?GO\:[0-9]{7,7}\)/g){
			my $go = $&;
			my $goid = substr($go,-11,-1);		# GO:####### を取り出す。
			if(exists($goterm_bioproc{$goid})){
				my $index = $goindex_bioproc{$goid};	# hash は遅そうなので取り出すのは１回にす	る
				$goenrich1_bioproc[$index] *= $ave_fc1[$g];
				$goenrich2_bioproc[$index] *= $ave_fc2[$g];			# foldchange 値は掛け算で積算する
				$goenrich3_bioproc[$index] *= $ave_fc3[$g];
				$goenrich4_bioproc[$index] *= $ave_fc4[$g];
				$gogenes_bioproc[$index]++;
			}else{
				$goindex_bioproc{$goid} = $gos_bioproc;			# 配列のindexを登録
				$goid_bioproc[$gos_bioproc] = $goid;			# GO番号 を登録
				$goterm_bioproc{$goid} = substr($go,20,-13);	# GO term を取り出す。
				$gogenes_bioproc[$gos_bioproc]=1;
				$goenrich1_bioproc[$gos_bioproc] = $ave_fc1[$g];
				$goenrich2_bioproc[$gos_bioproc] = $ave_fc2[$g];	# foldchange 値を４実験区分
				$goenrich3_bioproc[$gos_bioproc] = $ave_fc3[$g];
				$goenrich4_bioproc[$gos_bioproc] = $ave_fc4[$g];
				$gos_bioproc++;
			}
		}
	}
}

print 'Genes with GO term(s) in RAP3 annotation: '.$gene_with_gos."\n";
print 'GO tems found for Molecular Function: '.$gos_molfunc."\n";
print 'GO tems found for Cellular Component: '.$gos_celcomp."\n";
print 'GO tems found for Biological Process: '.$gos_bioproc."\n";

# 乗算したので遺伝子数分平方根する。ついでに２つの再現実験区の掛け算平方根も計算しておく
my @spkbc_molfunc;
my @camvc_molfunc;
my @spkbc_celcomp;
my @camvc_celcomp;
my @spkbc_bioproc;
my @camvc_bioproc;
for(my $i=0;$i<$gos_molfunc;$i++){
	$goenrich1_molfunc[$i] = $goenrich1_molfunc[$i]**(1/$gogenes_molfunc[$i]);
	$goenrich2_molfunc[$i] = $goenrich2_molfunc[$i]**(1/$gogenes_molfunc[$i]);
	$goenrich3_molfunc[$i] = $goenrich3_molfunc[$i]**(1/$gogenes_molfunc[$i]);
	$goenrich4_molfunc[$i] = $goenrich4_molfunc[$i]**(1/$gogenes_molfunc[$i]);
	$spkbc_molfunc[$i] = sqrt($goenrich1_molfunc[$i] * $goenrich2_molfunc[$i]);
	$camvc_molfunc[$i] = sqrt($goenrich3_molfunc[$i] * $goenrich4_molfunc[$i]);
}
for(my $i=0;$i<$gos_celcomp;$i++){
	$goenrich1_celcomp[$i] = $goenrich1_celcomp[$i]**(1/$gogenes_celcomp[$i]);
	$goenrich2_celcomp[$i] = $goenrich2_celcomp[$i]**(1/$gogenes_celcomp[$i]);
	$goenrich3_celcomp[$i] = $goenrich3_celcomp[$i]**(1/$gogenes_celcomp[$i]);
	$goenrich4_celcomp[$i] = $goenrich4_celcomp[$i]**(1/$gogenes_celcomp[$i]);
	$spkbc_celcomp[$i] = sqrt($goenrich1_celcomp[$i] * $goenrich2_celcomp[$i]);
	$camvc_celcomp[$i] = sqrt($goenrich3_celcomp[$i] * $goenrich4_celcomp[$i]);
}
for(my $i=0;$i<$gos_bioproc;$i++){
	$goenrich1_bioproc[$i] = $goenrich1_bioproc[$i]**(1/$gogenes_bioproc[$i]);
	$goenrich2_bioproc[$i] = $goenrich2_bioproc[$i]**(1/$gogenes_bioproc[$i]);
	$goenrich3_bioproc[$i] = $goenrich3_bioproc[$i]**(1/$gogenes_bioproc[$i]);
	$goenrich4_bioproc[$i] = $goenrich4_bioproc[$i]**(1/$gogenes_bioproc[$i]);
	$spkbc_bioproc[$i] = sqrt($goenrich1_bioproc[$i] * $goenrich2_bioproc[$i]);
	$camvc_bioproc[$i] = sqrt($goenrich3_bioproc[$i] * $goenrich4_bioproc[$i]);
}

open OUT, '>GOenrich_molfunc.txt';
print OUT 'GO ID'."\t";
print OUT 'GO term'."\t";
print OUT 'genes'."\t";
print OUT 'SpkBC1'."\t";
print OUT 'SpkBC2'."\t";
print OUT '35S-C1'."\t";
print OUT '35S-C2'."\t";
print OUT 'SpkBC'."\t";
print OUT '35S-C'."\n";
for(my $i=0;$i<$gos_molfunc;$i++){
	print OUT $goid_molfunc[$i]."\t";
	print OUT $goterm_molfunc{$goid_molfunc[$i]}."\t";
	print OUT $gogenes_molfunc[$i]."\t";
	print OUT $goenrich1_molfunc[$i]."\t";
	print OUT $goenrich2_molfunc[$i]."\t";
	print OUT $goenrich3_molfunc[$i]."\t";
	print OUT $goenrich4_molfunc[$i]."\t";
	print OUT $spkbc_molfunc[$i]."\t";
	print OUT $camvc_molfunc[$i]."\n";
}
close OUT;

open OUT, '>GOenrich_celcomp.txt';
print OUT 'GO ID'."\t";
print OUT 'GO term'."\t";
print OUT 'genes'."\t";
print OUT 'SpkBC1'."\t";
print OUT 'SpkBC2'."\t";
print OUT '35S-C1'."\t";
print OUT '35S-C2'."\t";
print OUT 'SpkBC'."\t";
print OUT '35S-C'."\n";
for(my $i=0;$i<$gos_celcomp;$i++){
	print OUT $goid_celcomp[$i]."\t";
	print OUT $goterm_celcomp{$goid_celcomp[$i]}."\t";
	print OUT $gogenes_celcomp[$i]."\t";
	print OUT $goenrich1_celcomp[$i]."\t";
	print OUT $goenrich2_celcomp[$i]."\t";
	print OUT $goenrich3_celcomp[$i]."\t";
	print OUT $goenrich4_celcomp[$i]."\t";
	print OUT $spkbc_celcomp[$i]."\t";
	print OUT $camvc_celcomp[$i]."\n";
}
close OUT;

open OUT, '>GOenrich_bioproc.txt';
print OUT 'GO ID'."\t";
print OUT 'GO term'."\t";
print OUT 'genes'."\t";
print OUT 'SpkBC1'."\t";
print OUT 'SpkBC2'."\t";
print OUT '35S-C1'."\t";
print OUT '35S-C2'."\t";
print OUT 'SpkBC'."\t";
print OUT '35S-C'."\n";
for(my $i=0;$i<$gos_bioproc;$i++){
	print OUT $goid_bioproc[$i]."\t";
	print OUT $goterm_bioproc{$goid_bioproc[$i]}."\t";
	print OUT $gogenes_bioproc[$i]."\t";
	print OUT $goenrich1_bioproc[$i]."\t";
	print OUT $goenrich2_bioproc[$i]."\t";
	print OUT $goenrich3_bioproc[$i]."\t";
	print OUT $goenrich4_bioproc[$i]."\t";
	print OUT $spkbc_bioproc[$i]."\t";
	print OUT $camvc_bioproc[$i]."\n";
}
close OUT;

# t値からP値を計算する関数。自由度６のt分布を端まで積分する。
sub p_value{
	my $argv1 = $_[0];
	my $start = $argv1;

	my $result=0;
	for(my $i=$start;$i<$t_end;$i+=$t_step){
		my $left = $t_dist_factor*(1+($i**2)/6)**(-3.5);
		my $right = $t_dist_factor*(1+(($i+$t_step)**2)/6)**(-3.5);
		my $area = ($left+$right)*$t_step/2;
		$result += $area;
	}
	$result = $result/$t_area_total;
	return($result);
}

