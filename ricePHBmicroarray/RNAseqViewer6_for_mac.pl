##################
#
# RNAseq Viewer version 2019.4.23 by 草野博彰＠京都大学生存圏研究所 【部外秘・無断配布禁止】
#
##################
#
# RNAseq Viewer は de Novo RNAseq 解析 で得られたデータを解析するために製作されました。
# Trinity + RSEM で解析したデータをBLAST やキーワードで検索し、散布図と画像で表現します。
# NCBI-SRA 等からダウンロードできるような生リードデータの解析はできません。
#
# BLASTかキーワードでコンティグを探し発現量比の分布を示す散布図にプロットすることができます。
# 散布図上のプロットをクリックすることでもコンティグの検索が可能です。
# 検索したコンティグを元に再度BLAST検索することで類似コンティグを探したりするために使います。
#
# 【 キーワード検索 】
# Keyword Search ボタンの横にキーワードを入力してKeyword Search ボタンを押す。
# RSEMデータに含まれる情報からキーワードにヒットするcontigを抜き出してscatter plot 表示します。
#
# 【 Local Blast 検索 】
# ローカルライブラリを検索します。alignment 、scatter plot、塩基配列リストに反映されます。
# 準備：ローカルライブラリが必要です。左下のPrepare Local BLAST library で検索対象にするデータベースの元データ（fasta形式）を選択します。
# 準備：ローカルライブラリは一度作れば二回作る必要はありません。左上の大きいLoadボタンを押すと記憶されます。
# load query ボタンを押してFasta形式のDNAまたはタンパク質の配列データをロードしLocal BLAST ボタンを押す
# または、load query の上の窓に塩基配列またはアミノ酸配列をコピー＆ペースト（ペーストはF4）してLocal BLAST ボタンを押す
# queryの種別(DNA/Protein) は自動認識されません。左下のチェックボックスで blastn / tblastn を選択する。
#
# 【 Tips 】
# コピーは command+C、ペーストは F4。command + V がなぜか機能しない。
# scatter plot の縦軸と横軸は Add New Plot ボタンの周りで設定します。
# 拡大縮小はそれぞれPlot、 Align ボタンの下のリストから拡大率を選択します。ボタンを押すと描画します
# 全プロット背景を含むscatter plot は画像ファイルに保存されます。ウィンドウでは全プロット背景が表示されません。
# scatter plot でクリックしたcontig の塩基配列は右上の窓に表示されます。この窓にcontig名を直接入力するとscatter plot に反映されます。
# 右下の窓はrevcom機能を備えた配列エディタです。この窓だけ折り返し表示されます。
#
# 【 システム要件 】
# ・macOSX 14 （本バージョンはmac用です）
# ・Perl 5.28.1に Tk、GD（pkgconfig, libgd依存）、Bioperl、LWP（HTML::Target依存）がインストールされていること
# ・NCBI blast+、XQuartz、Xcodes (runtime environment) がインストールされていること
# ・TrinityデータからBLASTn用のデータベースが構築されていること（ncbi blast+ に付いてくるmakeblastdb.exe を使う）
# 【 動作確認済みシステム 】
# ・macOSX mojave (10.14) 、Macbook Air および iMac （一体型デスクトップ）
#
##################

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Tools::Run::RemoteBlast;
use Tk;
use GD;
use MIME::Base64;		# GDで作った画像をtkに渡すために必要。なお、拡大縮小は失敗した
use Tk::PNG;			# 同上
my $ttfont = '/Library/Fonts/Arial.ttf';

my $workidentifier = int(rand(100000));

# メインの重たいデータファイルのパス
my $rsem_file = 'PR0614_result_ag.txt';		# RSEM データのファイル
my $rsem_file_status = 'not available';		# RSEM データの読み込み状況
my $trinity_file = 'Trinity.fasta';			# Trinity データのファイル
my $trinity_file_status = 'not available';	# Trinity データの読み込み状況

# blast結果とユーザー指定のfastaファイル
my $blast_file = 'localblast.txt';			# BLAST データのファイル
my $blast_file_status = 'not available';	# BLAST データの読み込み状況
my $entrylist_file = 'contiglist.fasta';		# IDリストのファイル
my $entrylist_file_status = 'not available';# IDリストの読み込み状況
my $ncbiblast_file = 'ncbiblast.txt';			# BLAST データのファイル
my $ncbiblast_file_status = 'not available';	# BLAST データの読み込み状況

# RSEM データを持つための変数
my @contig_fpkm;
my @contig_name;
my @contig_length;
my @contig_line;
my $contigs=0;
my $fpkms=0;
my @sample_name;
my $firstlinersem;

# Trinity データを持つための変数
my @trinity_seq;
my @trinity_name;
my $trinities=0;

# NCBI blast からe-valの積算値分布を読み出したデータをクロマトとして持つための変数
my @ncbi_chromat;
my $ncbi_chromat_max = 1;

# Local Blast から読み出したデータを持つための変数
my $query_length=0;	# Local blast に使ったquery の長さ。Local blast データがロードされたフラグとしても使う
my @hit_name;		# ヒットしたエントリのなまえ
my @subject_end;	# ヒットしたエントリーの長さ
my @subject_strand;	#
my $hits;			# ヒットしたエントリの数
my $ranks_max=50;	# アラインメントに表示するエントリの数
my $ranks=$ranks_max;	# アラインメントに表示するエントリの数
my @hsp_stt;		# セグメントの開始点	[エントリ][セグメント]	query側から見た数字
my @hsp_end;		# セグメントの終了点	[エントリ][セグメント]
my @hsp_stt_hit;	# セグメントの開始点	[エントリ][セグメント]	subject 側からみた数字
my @hsp_end_hit;	# セグメントの終了点	[エントリ][セグメント]
my @hsps;			# セグメント数			[エントリ]
my @cons_stt;	# 各セグメント内の連続一致領域の開始点	[エントリ][セグメント][連続一致部位の番号]
my @cons_end;	# 各セグメント内の連続一致領域の終了点	[エントリ][セグメント][連続一致部位の番号]
my @cons;		# 各セグメント内の連続一致領域の数		[エントリ][セグメント]
my @ord;
my @ords;
my $align_zoom = 0.5;
my $align_savefile;	# 画像をユーザー指定のファイルに保存するためのファイル名。デフォ値はundef。画像は自動的にテンポラリファイルに保存される。
		$align_savefile ='RNAseqViewer_blastalign.png';
my $switch_local=0;	# local blast を結果のファイルから呼ぶか今から実行するかのスイッチ

# 散布図画像を持つための変数
my $reset_scat_back=1;		# 散布図を描き直す必要が生じたとき
my $backimage;				# イメージオブジェクト。
my $saveimage_scat;
my $graphsize=900;			# キャンバスいっこ分のサイズ。拡縮するとき動的にするか？
my $scat_savefile;			# 画像をユーザー指定のファイルに保存するためのファイル名。デフォ値はundef。画像は自動的にテンポラリファイルに保存される。
		$scat_savefile ='RNAseqViewer_scatplot.png';
my $scat_showtext_limit = 5;

# ユーザー指定のリスト
my @fasta_name;
my @fasta_seq;
my @fasta_fpkm;
my $fastas=0;


# 散布図を Local blast でハイライトするための変数
my @hitcontig_fpkm;		# FPKM の値はRSEMから。
my $hitcontig_resetflag = 1;
my @hitcontig_name;		# 他はLocalblastからロードする。
my @hitcontig_length;
my $hitcontig_length_max;
my $hitcontig_length_min;
my @hitcontig_score;
my $hitcontig_score_max;
my $hitcontig_score_min;
my @hitcontig_evalue;
my @hitcontig_seq;
my $hitcontigs=0;
my @hitcontig_plotx;	# プロット位置。canvas上の座標。[][] 二番目の引数がプロットの番号
my @hitcontig_ploty;

# スキャッタープロットの数とXY軸と識別ラベル
my @scat_label;
$scat_label[0] = 'volcano';
$scat_label[1] = 'expr/cont';
my @scat_xfpkm;
my @scat_yfpkm;
$scat_xfpkm[0] = 2;$scat_yfpkm[0] = 3;
$scat_xfpkm[1] = 1;$scat_yfpkm[1] = 0;
$scat_xfpkm[2] = 4;$scat_yfpkm[2] = 5;
$scat_xfpkm[3] = 1;$scat_yfpkm[3] = 4;
$scat_xfpkm[4] = 1;$scat_yfpkm[4] = 5;
my $scats=5;
#my @scat_label;$scat_label[0] = 'expr/cont';$scat_label[1] = 'error';
#my @scat_xfpkm;$scat_xfpkm[0] = 0;$scat_xfpkm[1] = 5;
#my @scat_yfpkm;$scat_yfpkm[0] = 1;$scat_yfpkm[1] = 4;
#my $scats=2;
my $scat_newlabel;		# スキャッタープロットの入力用バッファ
my $scat_newx;
my $scat_newy;
my $scat_newx_buf = 'X'.$scat_newx;
my $scat_newy_buf = 'Y'.$scat_newy;
my $scat_hitextswitch = 1;
my $scat_plotcolor_switch = 0;
my $scat_zoom = 0.3;

# テキスト窓の変数
my $textwidget_status = 'empty';
my $textwidget_refreshflag = 0;

# Blast に投げる関連の変数
my $blastquery_status = 'invalid';
my $blast_query;
my $blastfactory;	# リモートblastのblastfactory
my $ncbiblast_status = 'RemoteBlast: idle';	# リモートblastの状態メッセージ
my $ncbiblast_status_count = 0;
my $local_library = 'Taxus180803';
my $local_library_status = 'not available';
my $blastmode_switch = 0;			# blastn か tblastn かどっちを使うかのスイッチ
my $blast_eval = '1e-2';			# local blast に与えるcutoff e-value

# キーワード検索関連の変数
my $keyword_input;
my $keyword_error = 'empty';
my @keyword_hit_name;
my @keyword_hit_seq;
my @keyword_hit_fpkm;
my @keyword_hit_length;
my $keyword_hits = 0;
my $keyword_hit_length_max;
my $keyword_hit_length_min;
my @keyword_hit_plotx;	# プロット位置。canvas上の座標。[][] １番目の引数がグラフの番号
my @keyword_hit_ploty;

# 前回のデータベースファイルの場所とかを思い出す
if( -f 'RNAseqViewer_settings.txt'){
	open IN, 'RNAseqViewer_settings.txt';		# 前回のセッティング情報を覚えておくようにセット
		$rsem_file = <IN>;chomp($rsem_file);
		$trinity_file =<IN>;chomp($trinity_file);
		$local_library = <IN>;chomp($local_library);
		if(-f $rsem_file){$rsem_file_status = 'available';}else{$rsem_file_status = 'not available';}
		if(-f $trinity_file){$trinity_file_status = 'available';}else{$trinity_file_status = 'not available';}
		$local_library_status = 'not available';
		if(-f $local_library.'.nhr'){$local_library_status = 'available';}
		if(-f $local_library.'.nin'){$local_library_status = 'available';}
		if(-f $local_library.'.nsq'){$local_library_status = 'available';}
		$scats = 0;
		while(my $line = <IN>){
			if($line =~ /[0-9]/){
				chomp $line;
				my @s = split(/\t/,$line);
				$scat_xfpkm[$scats] = $s[0];			# fpkmsより大きい数字だったときの対応はRSEMデータロード時にやるのでここではしない
				$scat_yfpkm[$scats] = $s[1];
				$scat_label[$scats] = $s[2];
				$scats++;
			}
		}
	close IN;
}

#############
# GUIを作る
#############

# ウィンドウを３つ
my $control_window = MainWindow -> new();
$control_window -> title('RNAseq Viewer '.$workidentifier);

my $alignment_window = $control_window -> Toplevel();
   $alignment_window -> title('Alignment View '.$workidentifier);
   $alignment_window -> geometry("700x300");
my $align_frame = $alignment_window -> Scrolled('Canvas',-scrollbars => 'se') -> pack(-fill => 'both', -expand => 1);
my $align_canvas = $align_frame -> Subwidget('scrolled') -> pack(-fill => 'both', -expand => 1);

my $plot_window = $control_window -> Toplevel();
   $plot_window -> title('ScatterPlot View '.$workidentifier);
   $plot_window -> geometry("400x400");
my $plot_frame = $plot_window -> Scrolled('Canvas',-scrollbars => 'se') -> pack(-fill => 'both', -expand => 1);
my $plot_canvas = $plot_frame -> Subwidget('canvas') -> pack(-fill => 'both', -expand => 1);


# 重たいデータをロードするボタンは上の方に置いておく
my $menuline1 = $control_window -> Frame() -> pack(-side => 'left', -anchor => 'nw');
	my $menuline11 = $menuline1 -> Frame() -> pack(-side => 'top', -anchor => 'nw');
		my $menuline112 = $menuline11 -> Frame() -> pack(-side => 'left');
			my $button_loaddata = $menuline112 -> Button(-text => "Load RNAseq Data\nand BLAST library\n\* wait a minute \*", -height => 5, -command => [\&load_data]) -> pack(-side => 'top');
		my $menuline111 = $menuline11 -> Frame() -> pack(-side => 'left');
			my $menuline1111 = $menuline111 -> Frame()-> pack(-anchor => 'nw');
				my $button_rsemfile = $menuline1111 -> Button(-text => 'RSEM data', -command => [\&findfile_rsem]) -> pack(-side => 'left');
				my $filepass_rsem   = $menuline1111 -> Entry(-textvariable => \$rsem_file, -width => 30) ->pack(-side => 'left');
				my $rsem_file_status_label = $menuline1111 -> Label(-textvariable => \$rsem_file_status) ->pack(-side => 'left');
#				my $load_rsem_button = $menuline1111 -> Button(-text => 'Load RSEM', -command => [\&load_rsem]) -> pack(-side => 'left');

			my $menuline1112 = $menuline111 -> Frame() -> pack(-anchor => 'nw');
				my $button_trinityfile = $menuline1112 -> Button(-text => 'Trinity data ', -command => [\&findfile_trinity]) -> pack(-side => 'left');
				my $filepass_trinity   = $menuline1112 -> Entry(-textvariable => \$trinity_file, -width => 30) ->pack(-side => 'left');
				my $trinity_file_status_label = $menuline1112 -> Label(-textvariable => \$trinity_file_status) ->pack(-side => 'left');
#				my $load_trinity_button = $menuline1112 -> Button(-text => 'Load Trinity ', -command => [\&load_trinity]) -> pack(-side => 'left');

			my $menuline1113 = $menuline111 -> Frame() -> pack(-anchor => 'nw');
				my $button_libraryfile = $menuline1113 -> Button(-text => ' BLAST lib  ', -command => [\&findfile_blastlibrary]) -> pack(-side => 'left');
				my $filepass_library   = $menuline1113 -> Entry(-textvariable => \$local_library, -width => 30) ->pack(-side => 'left');
				my $library_file_status_label = $menuline1113 -> Label(-textvariable => \$local_library_status) ->pack(-side => 'left');

	my $menuline12 = $menuline1 -> Frame() -> pack(-side => 'top', -anchor => 'nw');
		my $menuline121 = $menuline12 -> Frame() -> pack(-side => 'left', -anchor => 'nw');
			my $menuline1211 = $menuline121 -> Frame() -> pack(-side => 'left', -anchor => 'n');
				my $scatx_listbox_label = $menuline1211 -> Label(-textvariable => \$scat_newx_buf) ->pack(-side => 'top');
				my $scatx_listbox = $menuline1211 -> Scrolled('Listbox', -scrollbars => 'e', -height => 9, -width => 4)->pack(-side => 'left');
				$scatx_listbox -> bind("<Button-1>",\&select_newscatx);
			my $menuline1212 = $menuline121 -> Frame() -> pack(-side => 'left', -anchor => 'n');
				my $scaty_listbox_label = $menuline1212 -> Label(-textvariable => \$scat_newy_buf) ->pack(-side => 'top');
				my $scaty_listbox = $menuline1212 -> Scrolled('Listbox', -scrollbars => 'e', -height => 9, -width => 4)->pack(-side => 'left');
				$scaty_listbox -> bind("<Button-1>",\&select_newscaty);
				# とりあえず初期値を与えておく。あとで消すか
				for(my $i=0;$i<$scats;$i++){$scatx_listbox -> insert($i,$i);$scaty_listbox -> insert($i,$i);}
			my $menuline1213 = $menuline121 -> Frame() -> pack(-side => 'top', -anchor => 'nw');
				$menuline1213->Label(-text => 'Graph Title') -> pack(-side => 'top', -anchor => 'nw');
				my $scat_newlabel_entry = $menuline1213 -> Entry(-textvariable => \$scat_newlabel, -width => 15)->pack(-side => 'top', -anchor => 'n');
				my $scat_addbutton = $menuline1213 -> Button(-text => 'Add New Plot', -command => [\&add_new_plot], -height => 2) ->pack(-side => 'top', -anchor => 'ne', -fill => 'x');
				$menuline1213 -> Label(-text => 'Plot Coloring') -> pack(-side => 'top', -anchor => 'nw');
				$menuline1213 -> Radiobutton(-text => 'Blast Score',   -variable => \$scat_plotcolor_switch, -value => 0)-> pack(-side => 'top', -anchor => 'nw');
				$menuline1213 -> Radiobutton(-text => 'Blast e-value', -variable => \$scat_plotcolor_switch, -value => 1)-> pack(-side => 'top', -anchor => 'nw');
				$menuline1213 -> Radiobutton(-text => 'Contig length', -variable => \$scat_plotcolor_switch, -value => 2)-> pack(-side => 'top', -anchor => 'nw');
		my $scat_listbox = $menuline12 -> Scrolled( 'Listbox', -scrollbars => 'se', -height => 9, -width =>20 )->pack(-side => 'left');
			for(my $i=0;$i<$scats;$i++){$scat_listbox -> insert($i, 'x#'.$scat_xfpkm[$i].', y#'.$scat_yfpkm[$i].', '.$scat_label[$i]);}
			$scat_listbox -> bind("<Button-1>",\&select_scat);	# 左クリック
			$scat_listbox -> bind("<Button-3>",\&delete_scat);	# 右クリック

		my $menuline122 = $menuline12 -> Frame()-> pack(-side => 'left', -anchor => 'nw');
			my $menuline1224 = $menuline122 -> Frame()-> pack(-side => 'left', -anchor => 'nw');
				my $button_draw = $menuline1224 -> Button(-text => "Plot", -command => [\&draw_scat], -height => 1, -width => 4) -> pack(-side => 'top');
				$menuline1224 -> Label(-textvariable => \$scat_zoom)->pack(-side => 'top');
				my $zoomscat_listbox = $menuline1224 -> Scrolled('Listbox', -scrollbars => 'e', -width =>4, -height =>7) -> pack(-side => 'top');
				my @zoomscat_defaultlist = (0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0);
				for(my $i=1;$i<9;$i++){$zoomscat_listbox -> insert($i,$zoomscat_defaultlist[$i-1]);}
				$zoomscat_listbox -> bind("<Button-1>",\&set_scat_zoom);
			my $menuline1223 = $menuline122 -> Frame() -> pack(-side => 'left', -anchor => 'nw');
				$menuline1223 -> Button(-text => "Align", -height => 1, -width =>4, -command => [\&draw_align]) -> pack(-side => 'top');
				$menuline1223 -> Label(-textvariable => \$align_zoom)->pack(-side => 'top');
				my $zoomalign_listbox = $menuline1223 -> Scrolled('Listbox', -scrollbars => 'e', -width =>4, -height =>7) -> pack(-side => 'top');
				my @zoomalign_defaultlist = (0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0);
				for(my $i=1;$i<8;$i++){$zoomalign_listbox -> insert($i,$zoomalign_defaultlist[$i-1]);}
				$zoomalign_listbox -> bind("<Button-1>",\&set_align_zoom);

	my $menuline41 = $menuline1 -> Frame() -> pack(-anchor => 'nw', -side => 'left');
		my $menuline411 = $menuline41 -> Frame() -> pack(-anchor => 'nw');
			my $button_blastfile = $menuline411 -> Button(-text => 'Load Local BLAST result file', -command => [\&findfile_blast]) -> pack(-side => 'left');
			my $filepass_blast   = $menuline411 -> Entry(-textvariable => \$blast_file, -width => 16) ->pack(-side => 'left');
			my $blast_file_status_label = $menuline411 -> Label(-textvariable => \$blast_file_status) ->pack(-side => 'left');
			my $button_localblastfile_clear = $menuline411 -> Button(-text => 'clear', -command => [\&clear_localblast]) -> pack(-side => 'left');

		my $menuline414 = $menuline41 -> Frame() -> pack(-anchor => 'nw');
			my $button_ncbiblastfile = $menuline414 -> Button(-text => 'Load NCBI BLAST result file ', -command => [\&findfile_ncbiblast]) -> pack(-side => 'left', -anchor => 'nw');
			my $filepass_ncbiblast   = $menuline414 -> Entry(-textvariable => \$ncbiblast_file, -width => 16) ->pack(-side => 'left');
			my $ncbiblast_file_status_label = $menuline414 -> Label(-textvariable => \$ncbiblast_file_status) ->pack(-side => 'left');
			my $button_ncbiblastfile_clear = $menuline414 -> Button(-text => 'clear', -command => [\&clear_ncbiblast]) -> pack(-side => 'left');
		my $menuline415 = $menuline41 -> Frame() -> pack(-anchor => 'nw');
			my $button_createdb = $menuline415 -> Button(-text => 'Prepare Local BLAST library ', -command => [\&create_localblastdb]) -> pack(-side => 'left', -anchor => 'nw');
			my $filepass_createdb   = $menuline415 -> Entry(-textvariable => \$local_library, -width => 16) ->pack(-side => 'left');
			my $createdb_label = $menuline415 -> Label(-textvariable => \$local_library_status) ->pack(-side => 'left');
		my $menuline416 = $menuline41 -> Frame() -> pack(-anchor => 'nw');
				$menuline416 -> Label(-text => 'select a Local Blast Mode') -> pack(-side => 'top', -anchor => 'nw');
				$menuline416 -> Radiobutton(-text => 'blastn: use DNA seq query for Local BLAST',   -variable => \$blastmode_switch, -value => 0)-> pack(-side => 'top', -anchor => 'nw');
				$menuline416 -> Radiobutton(-text => 'tblastn: use Protein seq query for Local BLAST', -variable => \$blastmode_switch, -value => 1)-> pack(-side => 'top', -anchor => 'nw');
	my $menuline42 =$menuline1 -> Frame() -> pack(-side => 'left');
		$menuline42 -> Label(-textvariable => \$blast_eval)->pack(-side => 'top');
		my $seteval_listbox = $menuline42 -> Scrolled('Listbox', -scrollbars => 'e', -width =>4, -height =>7) -> pack(-side => 'top');
		my @seteval_defaultlist = (1e-2, 1e-5, 1e-10, 1e-20, 1e-40, 1e-60, 1e-180);
		for(my $i=1;$i<8;$i++){$seteval_listbox -> insert($i,$seteval_defaultlist[$i-1]);}
		$seteval_listbox -> bind("<Button-1>",\&set_blast_eval);


my $menuline2 = $control_window -> Frame() -> pack(-anchor => 'nw', -side => 'left');
	my $menuline21 = $menuline2 -> Frame() -> pack(-side => 'top', -anchor => 'nw', -fill => 'x');
		my $menuline211 = $menuline21 -> Frame() -> pack(-side => 'left', -anchor => 'w');
			my $keyword_entrybox   = $menuline211 -> Entry(-textvariable => \$keyword_input, -width => 25) ->pack(-side => 'top');
			   $keyword_entrybox -> bind("<KeyRelease>", \&keyword_errorcheck);
			my $button_keyword_clear = $menuline211 -> Button(-text => 'clear', -command => [\&clear_keyword]) -> pack(-side => 'right');
			my $error_keyword_search = $menuline211 -> Label(-textvariable => \$keyword_error) ->pack(-anchor => 'w', -side => 'left');
		my $button_keyword_search = $menuline21 -> Button(-text => "Keyword\nSearch", -height => 3, -width => 10, -command => [\&search_keyword]) -> pack(-anchor => 'nw', -side => 'left');
		my $button_exlocalblast  = $menuline21 -> Button(-text => "Local\nBLAST", -height => 3, -width => 10, -command => [\&exec_localblast]) -> pack(-side => 'right', -anchor => 'nw');
		my $menuline212 = $menuline21 -> Frame() -> pack(-side => 'right', -anchor => 'se');
			my $button_exremote_submit = $menuline212 -> Button(-text => " Submit", -height => 1, -command => [\&exec_remoteblast]) -> pack(-side => 'top', -anchor => 'nw');
			my $button_exremote_retrieve = $menuline212 -> Button(-text => "retrieve", -height => 1, -command => [\&retrieve_remoteblast]) -> pack(-side => 'top', -anchor => 'nw');

	my $menuline22 = $menuline2 -> Frame() -> pack(-side => 'top', -anchor => 'nw');
		my $menuline221 = $menuline22 -> Frame() -> pack(-side => 'top', -anchor => 'nw');
			my $blast_query_input = $menuline221 -> Entry(-textvariable => \$blast_query, -width => 60) ->pack(-side => 'top');
				$blast_query_input -> bind("<KeyRelease>", \&blastquery_refresh);
				$blast_query_input -> eventAdd(qw[<<Paste>> <Key-F5>]);
			my $button_blastqueryfasta = $menuline221 -> Button(-text => 'load query', -command => [\&findfile_blastquery]) -> pack(-side => 'left');
			my $blastquery_file_status_label = $menuline221 -> Label(-textvariable => \$blastquery_status) ->pack(-anchor => 'e', -side => 'left');
			my $button_clearblastqueryfasta = $menuline221 -> Button(-text => 'clear', -command => [\&clear_blastquery]) -> pack(-anchor => 'w',-side => 'right');
			my $revcom_button = $menuline221 -> Button(-text => 'revcom', -command => [\&revcom_query])->pack(-side => 'right', -anchor => 'e');

#my $menuline3 = $menuline0 -> Frame() -> pack(-anchor => 'nw', -side => 'left');

my $menuline4 = $menuline2 -> Frame() -> pack(-anchor => 'nw', -side => 'left');
	my $textwidget2 = $menuline4 -> Scrolled('Text',-scrollbars => 'se', -font => 'Arial 12', -width => 45, -height => 16, -wrap => 'none' )->pack(-side => 'top', -fill => 'both');
#		my $menuline222 = $menuline4 -> Frame() -> pack(-side => 'top', -anchor => 'nw');
			my $button_localblastfile_clear = $menuline4 -> Button(-text => 'clear localblast', -command => [\&clear_localblast]) -> pack(-side => 'right');
			my $button_makelist = $menuline4 -> Button(-text => 'Search seq', -command => [\&blast2fasta]) -> pack(-side => 'left');

my $menuline5 = $control_window -> Frame() -> pack(-side => 'left', -anchor => 'nw');
	# blast結果からユーザー指定リストを作ってくれる便利機能を作ろう
	my $menuline53 = $menuline5 -> Frame()->pack(-anchor => 'nw');
	# ユーザー指定リストを表示する部分。Fastaをテキスト編集できるようにしたら便利か
	my $menuline55 = $menuline5 -> Frame()->pack(-side => 'top', -anchor => 'nw');
		my $highlightmessage1 = $menuline55 -> Label(-text => 'Show Contig Name:')-> pack(-side => 'left', -anchor => 'nw');
		my $highlight_on =  $menuline55 -> Radiobutton(-text => 'enable' , -variable => \$scat_hitextswitch, -value => 1)-> pack(-side => 'left', -anchor => 'nw');
		my $highlight_off = $menuline55 -> Radiobutton(-text => 'disable', -variable => \$scat_hitextswitch, -value => 0)-> pack(-side => 'left', -anchor => 'nw');
		my $highlightmessage2 = $menuline55 -> Label(-text => ' Limit:')-> pack(-side => 'left', -anchor => 'nw');
		my $highlightlimit_entry = $menuline55 -> Entry(-textvariable => \$scat_showtext_limit, -width => 3) ->pack(-side => 'left');

	my $textwidget = $menuline5 -> Scrolled('Text',-scrollbars => 'se', -font => 'Arial 12', -width => 40, -height => 13, -wrap => 'none' )->pack(-side => 'top', -fill => 'both');
	$textwidget -> bind("<KeyRelease>", \&textwidget_refresh);
	my $menuline6 = $menuline5->Frame()->pack(-side => 'top', -fill => 'x' );
		my $button_blastfile = $menuline6 -> Button(-text => 'Load', -command => [\&findfile_fasta]) -> pack(-side => 'left');
		my $savebutton_text = $menuline6 -> Button(-text => 'Save', -command => [\&save_fasta]) -> pack(-side => 'left');
		my $checkbutton_text = $menuline6 -> Button(-text => 'Search', -command => [\&input_fasta]) -> pack(-side => 'left');
		my $textwidget_status_display = $menuline6 -> Label(-textvariable => \$textwidget_status) -> pack(-side => 'left');
	my $clearbutton_text = $menuline6 -> Button(-text => 'clear', -command => [\&clear_fasta]) -> pack(-side => 'right');
	my $textwidget3 = $menuline5 -> Scrolled('Text',-scrollbars => 'e', -font => 'Arial 10', -width => 40, -height => 7, -wrap => 'char' )->pack(-side => 'top', -fill => 'both');
	my $button_ncbiblastfile_clear = $menuline5 -> Button(-text => 'clear chromat', -command => [\&clear_ncbiblast]) -> pack(-side => 'left');
	my $remoteblast_status_label = $menuline5 -> Label(-textvariable => \$ncbiblast_status) ->pack(-anchor => 'w', -side => 'left');
	my $revcom_button3 = $menuline5 -> Button(-text => 'revcom', -command => [\&revcom_free])->pack(-side => 'right');







MainLoop();
exit;

##################################################################################################### ここまでコントロールパネル  以下サブルーチン

#########
# Fastaファイルをテキストに表示するようにしたい。とりあえずサブ
#########
sub show_text(){
	my $t = $textwidget -> get('1.0','end-1c');
	print STDERR $t."\n";
}

#########
# Fastaファイルをテキストに表示する窓が変更されたときのフラグ
#########
sub textwidget_refresh(){
	$textwidget_refreshflag = 1;
	$textwidget_status = 'user edit';
}

#########
# 散布図設定のリスト画面が選択されたときの動作。特に何もしない。位置を黒窓に表示するだけ。
#########
sub select_scat(){
	my $t = $scat_listbox -> get('anchor');
	print STDERR 'selected: '.$t."\n";
}

#########
# アラインメント図のズーム値を設定するリストボックスの動作
#########
sub set_align_zoom(){
	my $t = $zoomalign_listbox -> get('anchor');
	$align_zoom = $t;
}

#########
# スキャッタープロット図のズーム値を設定するリストボックスの動作
#########
sub set_scat_zoom(){
	my $t = $zoomscat_listbox -> get('anchor');
	$scat_zoom = $t;
	$reset_scat_back = 1;
}

#########
# スキャッタープロットの設定で右クリックしたときの動作。該当リストを削除する
#########
sub delete_scat(){
	my $t = $scat_listbox -> get('anchor');
	if(defined($t)){
		my $user = $control_window -> messageBox (-type => 'okcancel', -message => 'Are you sure to exclude:'."\n".$t, -icon => 'info');
		if($user eq 'Ok'){
			my $n = $scat_listbox -> index('anchor');
			my @d = splice(@scat_label,$n,1);
			splice(@scat_xfpkm,$n,1);
			splice(@scat_yfpkm,$n,1);
			$scats--;
			foreach(@d){print 'delete '.$_."\n";}
			#foreach(@scat_label){print 'present '.$_."\n";}
			$scat_listbox -> delete('anchor');
			$reset_scat_back = 1;
		}
	}
}
#########
# 新しいスキャッタープロットの設定のxy値が選択されたときの動作 $scat_xxx[]で散布図の要素x,y,labelを持ち、$scatsで散布図の数を持つ。ユーザーへはリストボックスで表示
#########
sub select_newscatx(){
	$scat_newx = $scatx_listbox -> get('anchor');
	$scat_newx_buf = 'X='.$scat_newx;
	print STDERR $scat_newx_buf." selected\n";
}
sub select_newscaty(){
	$scat_newy = $scaty_listbox -> get('anchor');
	$scat_newy_buf = 'Y='.$scat_newy;
	print STDERR $scat_newy_buf." selected\n";
}
#########
# 新しいスキャッタープロットの設定を追加するボタンの動作
#########
sub add_new_plot(){
	if(defined($scat_newx)){
	if(defined($scat_newy)){
		my $check =1;
		for(my $i=0;$i<$scats;$i++){
			if($scat_xfpkm[$i] == $scat_newx){
			if($scat_yfpkm[$i] == $scat_newy){
			if($scat_label[$i] eq $scat_newlabel){
				$check = 0;
			}}}
		}
		if($check == 1){
			$scat_xfpkm[$scats] = $scat_newx;
			$scat_yfpkm[$scats] = $scat_newy;
			$scat_label[$scats] = $scat_newlabel;
			$scats++;
			$scat_listbox -> delete(0,'end');
			for(my $i=0;$i<$scats;$i++){$scat_listbox -> insert($i, 'x#'.$scat_xfpkm[$i].', y#'.$scat_yfpkm[$i].', '.$scat_label[$i]);}
			$reset_scat_back = 1;
		}
		for(my $i=0;$i<$fpkms;$i++){	# サンプル名を表示する
			print $i."\t".$sample_name[$i]."\n";
		}
		save_settings();
	}}
}

#########
# キーワード検索機能関連。local blast 関連の機能と表示窓など共有なので注意。
#########

sub clear_keyword(){
	$keyword_input = '';
	$keyword_error = 'empty';
	if($query_length == 0){
		# 検索結果の方のテキストを削除してすきゃったプロットするとか検索結果の配列を初期化するとか
		$textwidget2 -> delete('1.0','end');
		@keyword_hit_name = ();
		@keyword_hit_seq = ();
		@keyword_hit_fpkm = ();
		@keyword_hit_length = ();
		$keyword_hits = 0;

		draw_scat();
	}
}

#########
# キーワード検索。
#########
sub search_keyword(){
	if($query_length == 0){
	if($keyword_error eq 'OK'){
		@keyword_hit_name = ();
		@keyword_hit_seq = ();			# 塩基配列を探すのはボタンを押してからにすることにした
		@keyword_hit_fpkm = ();
		@keyword_hit_length = ();
		$keyword_hits = 0;
		
		my $filename_key = $keyword_input;
		$filename_key =~ s/\|/_/g;
		$filename_key =~ s/\ /_/g;
		$filename_key =~ s/\(/_/g;
		$filename_key =~ s/\)/_/g;
		$filename_key =~ s/\,/_/g;
		open OUT, '>temp_keyhits'.$workidentifier.'_'.$filename_key.'.txt';
		print OUT $firstlinersem;
		for(my $i=0;$i<$contigs;$i++){
			if($contig_line[$i] =~ /$keyword_input/){
				print OUT $contig_line[$i]."\n";
				$keyword_hit_name[$keyword_hits] = $contig_name[$i];
				$keyword_hit_length[$keyword_hits] = $contig_length[$i];
				for(my $j=0;$j<$fpkms;$j++){$keyword_hit_fpkm[$keyword_hits][$j] = $contig_fpkm[$i][$j];}
				$keyword_hits++;
			}
		}
		close OUT;
		print STDERR 'keyword search: '.$keyword_hits." contigs hit\n";
		$keyword_error = $keyword_hits.' contigs hit';		# ヒットしたコンティぐ件数を表示する。配列の検索はユーザーが判断する

		# length が長い順に並べる
		if($keyword_hits>0){
			my @temp_name;
			my @temp_length;
			my @temp_fpkm;
			for(my $i=0;$i<$keyword_hits;$i++){
				my $max = 0;
				my $maxid;
				for(my $j=0;$j<$keyword_hits;$j++){
					if($max < $keyword_hit_length[$j]){
						$max = $keyword_hit_length[$j];
						$maxid = $j;
					}
				}
				$temp_name[$i] = $keyword_hit_name[$maxid];
				$temp_length[$i] = $keyword_hit_length[$maxid];
				for(my $k=0;$k<$fpkms;$k++){$temp_fpkm[$i][$k] = $keyword_hit_fpkm[$maxid][$k];}
				$keyword_hit_length[$maxid] = 0;
			}
			for(my $i=0;$i<$keyword_hits;$i++){
				$keyword_hit_name[$i] = $temp_name[$i];
				$keyword_hit_length[$i] = $temp_length[$i];
				for(my $k=0;$k<$fpkms;$k++){$keyword_hit_fpkm[$i][$k] = $temp_fpkm[$i][$k];}
			}
			$keyword_hit_length_max = $keyword_hit_length[0];
			$keyword_hit_length_min = $keyword_hit_length[$keyword_hits-1];
		}

		# blast結果を表示する窓に表示する
		if($keyword_hits>0){
			$textwidget2 -> delete('1.0','end');
			for(my $i=0;$i<$keyword_hits;$i++){
				$textwidget2 -> insert('end','>'.$keyword_hit_name[$i]."\n");
				$textwidget2 -> insert('end',$keyword_hit_seq[$i]."\n");
			}
		}else{
			$textwidget2 -> delete('1.0','end');
			$textwidget2 -> insert('end','Not hits found for keyword: '.$keyword_input."\n");
		}
		draw_scat();
	}}
}

#########
# キーワード入力欄を触ったときの操作。キーをチェックする
#########
sub keyword_errorcheck(){

	if($query_length == 0){
		my $invalid_key = 0;
#		if($keyword_input =~ /\|/){$invalid_key = 1;}
		if($keyword_input =~ /\(/){$invalid_key = 1;}
		if($keyword_input =~ /\)/){$invalid_key = 1;}
#		if($keyword_input =~ /\,/){$invalid_key = 1;}
		if($keyword_input =~ /\&/){$invalid_key = 1;}
		if($keyword_input =~ /\"/){$invalid_key = 1;} #"
		if($keyword_input =~ /\!/){$invalid_key = 1;}
		if($keyword_input =~ /\\/){$invalid_key = 1;}

		if($invalid_key == 1){
			$keyword_error = 'invalid';
		}else{
			$keyword_error = 'OK';
			if($keyword_input eq ''){$keyword_error = 'empty';}
		}
	}else{
		$keyword_error = 'disable';
	}
}

#########
# Local BLAST のライブラリーを作る
#########

sub create_localblastdb(){
	my $buffer = $control_window -> getOpenFile( -filetypes => [['Fasta','.txt'],['Fasta','.fasta'],['Fasta','.fas'],['All files','*'] ]);
	if( -f $buffer){
		my $buffer_ou = substr($buffer, 0,length($buffer)-4);
		my $buffer_out = $buffer_ou.'.nin';
		if($buffer =~ /\.fasta/){
			$buffer_ou = substr($buffer, 0,length($buffer)-6);
			$buffer_out = $buffer_ou.'.nin';
		}
		if( -f $buffer_out){
			print STDERR 'the Local BLAST Library exists already: '.$buffer_ou."\n";
			$local_library = $buffer_ou;
			$local_library_status = 'available';
		}else{
			print STDERR 'create a Local BLAST Library: '.$buffer_ou."\n";
			system 'makeblastdb -in '.$buffer.' -dbtype nucl -out '.$buffer_ou;
			if( -f $buffer_out){
				$local_library = $buffer_ou;
				$local_library_status = 'available';
			}
		}
	}
}

#########
# Local BLAST のe-value値をセットする
#########

sub set_blast_eval(){
	my $t = $seteval_listbox -> get('anchor');
	$blast_eval = $t;
	$reset_scat_back = 1;
}

#########
# Local BLAST のローカルライブラリーへのパスを選ばせる
#########

sub findfile_blastlibrary(){
	my $buffer = $control_window -> getOpenFile( -filetypes => [['Local Library','.nhr'],['Local Library','.nin'],['Local Library','.nsq'],['All files','*'] ]);
	if( -f $buffer){
		$buffer = substr($buffer, 0,length($buffer)-4);
		print STDERR 'set Local BLAST Library: '.$buffer."\n";
		$local_library = $buffer;
		$local_library_status = 'available';
	}
}

#########
# BLAST に投げる配列をファイルからロードする操作
#########

sub findfile_blastquery(){
	my $buffer = $control_window -> getOpenFile( -filetypes => [['Fasta','.txt'],['Fasta','.fasta'],['Fasta','.fas'],['All files','*'] ]);
	if( -f $buffer){
		my $seqin = Bio::SeqIO->new(-file => $buffer, -format => 'fasta');
		if(my $seqobj = $seqin -> next_seq){
			$blast_query = $seqobj -> seq;
			$blastquery_status = 'available';
		}else{
			$blastquery_status = 'invalid';
		}
	}else{
		print STDERR 'canceled'."\n";
	}
}

#########
# 入力済みのBLAST クエリをクリアするボタンの動作
#########

sub clear_blastquery(){
	$blast_query = '';
	$blastquery_status = 'invalid';
}

#########
# BLAST クエリが手入力されたときの動作
#########

sub blastquery_refresh(){
	if($blast_query =~/\>/){
		$blastquery_status = 'invalid';	# fastaを貼った場合はどこまでが名前か判らないのでNGとする。
	}else{
		if($blast_query =~ /B|[D-F]|[H-M]|[O-S]|[U-Z]|b|[d-f]|[h-m]|[o-s]|[u-z]/){
			$blastquery_status = 'invalid';
		}else{
			$blast_query =~ s/\s//g;	# 改行とか空白とか入ってたら消す
			$blast_query =~ s/[0-9]//g;	# 数字も消す
			$blastquery_status = 'available';	# fasta の > が無ければいいということにしようか。本当はATGC以外の記号不可としたい。
	}}
}

#########
# リモートblastを実行する
#########

sub exec_remoteblast(){
	if(defined($blastfactory)){}else{
	if(length($blast_query)>0){
		open OUT,'>temp_blastquery_'.$workidentifier.'.txt';
		print OUT '>query'."\n".$blast_query;
		close OUT;
		@ncbi_chromat = ();
		$ncbi_chromat_max = 1;
		$blastfactory = new Bio::Tools::Run::RemoteBlast(-prog => 'blastx', -data => 'nr', -expect => '1e-10', -readmethod => 'Blast');
		print STDERR 'submit BLAST job to NCBI'."\n";
		$ncbiblast_status = 'RemoteBlast: blast job running';
		$ncbiblast_status_count = 1;
#		my $r = $blastfactory -> submit_blast('temp_blastquery_'.$workidentifier.'.txt');
		$blastfactory -> submit_blast('temp_blastquery_'.$workidentifier.'.txt');
	}}
}

#########
# リモートblastをretrieveする
#########
sub retrieve_remoteblast(){
	if(defined($blastfactory)){
#		BLA:while(my @rids = $blastfactory -> each_rid){
		my @rids = $blastfactory -> each_rid;
			foreach my $rid(@rids){
				my $rc = $blastfactory -> retrieve_blast($rid);
				if(!ref($rc)){
					if($rc<0){$blastfactory -> remove_rid($rid);}
					print STDERR 'Request to retrieve NCBI blast result remote: '.$ncbiblast_status_count."\n";
					$ncbiblast_status = 'RemoteBlast: retrieve request '.$ncbiblast_status_count.' times';
					$ncbiblast_status_count++;
				}else{
					print STDERR 'remoteblast done'."\n";
					$blastfactory -> remove_rid($rid);
					$blastfactory->save_output('temp_ncbiblast_'.$workidentifier.'.txt');
					$ncbiblast_file = 'temp_ncbiblast_'.$workidentifier.'.txt';
					$ncbiblast_file_status = 'available';
					my $blastresult_ncbi = $rc->next_result();

					@ncbi_chromat = ();
					$ncbi_chromat_max = 1;
					my $cutoff = 10;
					NCBI:while(my $hit = $blastresult_ncbi -> next_hit){
						my $judge = $hit -> significance;
						if($judge > 0){$judge = (-log($judge));}
						if($judge == 0){$judge = 181;}
						if($judge < $cutoff){last NCBI;}
						while (my $hsp = $hit -> next_hsp){
							my $start = $hsp -> start('query');
							my $end = $hsp -> end('query');
							my $eval = $hsp -> evalue;
							if($eval > 0){$eval = (-log($eval));}
							if($eval == 0){$eval = 181;}
							$ncbi_chromat[$start] += $eval;	if($ncbi_chromat_max < $ncbi_chromat[$start]){$ncbi_chromat_max = $ncbi_chromat[$start];}
							$ncbi_chromat[$end] += $eval;	if($ncbi_chromat_max < $ncbi_chromat[$end]  ){$ncbi_chromat_max = $ncbi_chromat[$end];  }
						}
					}
					print STDERR 'NCBI blast: maxpeak = '.$ncbi_chromat_max."\n";
					$ncbiblast_status = 'RemoteBlast: idle';
					undef($blastfactory);
					draw_align();
				}
			}
		#}
	}
}

#########
# NCBI blast の結果をリセットする
#########

sub clear_ncbiblast(){
	$ncbiblast_file = '';
	$ncbiblast_file_status = 'not available';
	$ncbi_chromat_max = 0;
	draw_align();
}

#########
# Local blast の結果をリセットする
#########

sub clear_localblast(){

	if($query_length > 0){
		$blast_file = '';
		$blast_file_status = 'not available';
		$query_length = 0;
		$hitcontigs = 0;
		$textwidget2 -> delete('1.0','end');
		keyword_errorcheck();
		draw_align();
		draw_scat();
	}
}

#########
# フリー窓のデータをrevcomする
#########

sub revcom_free(){
	my $text = $textwidget3 -> get('1.0','end-1c');
	if(defined($text)){
		if($text =~ /\>/){
			my @ss = split(/\>/,$text);
			$text = $ss[1];
			$text =~ s/\r//g;
			my @s = split(/\n/,$text);
			my $id = $s[0];
			$text =~ s/$id//;
			$text =~ s/\n//g;
			my $seq = $text;
			$id =~ s/\>//;
			my $for_obj = Bio::Seq -> new(-id => $id, -seq => $seq);
			my $rev_obj = $for_obj -> revcom;
			my $r = $rev_obj -> seq;
			$textwidget3 -> delete('1.0','end');
			$textwidget3 -> insert('end','>'.$id."_complement\n");
			$textwidget3 -> insert('end',$r."\n");
		}else{
			$text =~ s/\n//g;
			$text =~ s/\r//g;
			$text =~ s/\s//g;
			my $for_obj = Bio::Seq -> new(-id => 'test', -seq => $text);
			my $rev_obj = $for_obj -> revcom;
			my $r = $rev_obj -> seq;
			$textwidget3 -> delete('1.0','end');
			$textwidget3 -> insert('end',$r."\n");
		}
	}
}

#########
# ユーザー入力のblastクエリをrevcomする
#########

sub revcom_query(){
	if(defined($blast_query)){
		my $for_obj = Bio::Seq -> new(-id => 'test', -seq => $blast_query);
		my $rev_obj = $for_obj -> revcom;
		my $r = $rev_obj -> seq;
		$blast_query = $r;
	}
}


#########
# ローカルでblastを実行する。ライブラリの指定はどうしようかな。リモートと違って外部exeを呼ぶので面倒くさい
#########
sub exec_localblast(){

	if(length($blast_query)>0){
		unlink 'temp_localblast_'.$workidentifier.'.txt';
		open OUT,'>temp_blastquery_'.$workidentifier.'.txt';
		print OUT '>query'."\n".$blast_query;
		close OUT;
		my $executable = 'blastn';
		if($blastmode_switch == 1){$executable = 'tblastn';}
		system $executable.' -query temp_blastquery_'.$workidentifier.'.txt -db '.$local_library.' -evalue '.$blast_eval.' -out temp_localblast_'.$workidentifier.'.txt';
		my $check =0;
		while($check == 0){if( -f 'temp_localblast_'.$workidentifier.'.txt'){$check = 1;}else{print STDERR '.';sleep 5;}}
		if($check==1){
			print STDERR 'done'."\n";
			$blast_file = 'temp_localblast_'.$workidentifier.'.txt';
			$switch_local = 0;
			findfile_blast_body();
		}
	}
}


#########
# ファイルからLocal BLAST をロードする操作。メインデータなので大変
#########

sub findfile_blast(){
	$switch_local = 1;
	findfile_blast_body();
}

#########
# Local BLAST をロードする操作。メインデータなので大変
#########
sub findfile_blast_body(){

	my $buffer;
	if($switch_local == 1){$buffer = $control_window -> getOpenFile( -filetypes => [['BLAST result','.txt'],['All files','*'] ]);}
	if($switch_local == 0){$buffer = $blast_file;}
	if( -f $buffer){
		open IN, $buffer;my $test = <IN>;close IN;
		if($test =~ /BLAST/){
			$blast_file = $buffer;
			$blast_file_status = 'available';

			# データを変数に格納する
			my $searchio = new Bio::SearchIO (-format => 'blast', -file => $blast_file);
			my $blastresult = $searchio -> next_result;
			$query_length = $blastresult -> query_length;
			$hits = $blastresult -> num_hits;
			if($hits > $ranks_max){$ranks = $ranks_max;}
			if($hits < $ranks_max){$ranks = $hits;}

			for(my $i=0;$i<$ranks;$i++){
				my $hit = $blastresult -> next_hit;
				$subject_end[$i] = $hit -> length;
				$hit_name[$i] = $hit -> name;
				$hsps[$i]=0;
				while (my $hsp = $hit -> next_hsp){

					$hsp_stt[$i][$hsps[$i]] = ($hsp->start('query'));
					$hsp_end[$i][$hsps[$i]] = ($hsp->end('query'));
					$hsp_stt_hit[$i][$hsps[$i]] = ($hsp->start('hit'));
					$hsp_end_hit[$i][$hsps[$i]] = ($hsp->end('hit'));
					$subject_strand[$i] = $hsp->strand('hit');
					my @inds_query = $hsp -> seq_inds('query','conserved',1);
					# この関数は面倒。ヒット領域の数字を文字列として 256-364 みたいに返す。配列で。一塩基だけの場合は数字のみ
					$cons[$i][$hsps[$i]]=0;
					foreach my $in (@inds_query){
						my $st;
						my $en;
						if($in =~ /\-/){
							my @s = split(/\-/,$in);
							$st = $s[0];
							$en = $s[1];
						}else{
							$st = $in;
							$en = $in;
						}
						$cons_stt[$i][$hsps[$i]][$cons[$i][$hsps[$i]]] = $st;
						$cons_end[$i][$hsps[$i]][$cons[$i][$hsps[$i]]] = $en;
						$cons[$i][$hsps[$i]]++;
					}
					$hsps[$i]++;
				}
			}


			# 似てないけど配列は存在するところは枠内に含めるための面倒な計算
			for(my $i=0;$i<$ranks;$i++){
				# まずソートする
				my @order;
				my $orders=0;				# hsp番号を返す。
				while($orders<$hsps[$i]){
					my $min = $subject_end[$i]+2000;
					for(my $j=0;$j<$hsps[$i];$j++){
						my $check=1;
						for(my $k=0;$k<$orders;$k++){if($j == $order[$k]){$check=0;}}
						if($check==1){
							if($min > $hsp_stt[$i][$j]){
								$min = $hsp_stt[$i][$j];
								$order[$orders] = $j;
							}
						}
					}
					$orders++;
				}
				# ソートしたの便利なので取っとく
				for(my $j=0;$j<$orders;$j++){
					$ord[$i][$j] = $order[$j];
				}
				$ords[$i] = $orders;

				if($subject_strand[$i] == 1){
					# 一番左を足す
					if($hsp_stt_hit[$i][$order[0]] > 1){
						$hsp_stt[$i][$order[0]] = $hsp_stt[$i][$order[0]] - $hsp_stt_hit[$i][$order[0]];
					}
					# 一番右を足す
					if($hsp_end_hit[$i][$order[$orders-1]] < $subject_end[$i]){
						$hsp_end[$i][$order[$orders-1]] = $hsp_end[$i][$order[$orders-1]] + $subject_end[$i] - $hsp_end_hit[$i][$order[$orders-1]];
					}
					# ギャップを埋める
					for(my $j=0;$j<$orders-1;$j++){
						if($hsp_end_hit[$i][$order[$j]] < $hsp_stt_hit[$i][$order[$j+1]]){
							$hsp_end[$i][$order[$j]] = $hsp_end[$i][$order[$j]] + $hsp_stt_hit[$i][$order[$j+1]] - $hsp_end_hit[$i][$order[$j]];
						}
					}
				}else{
					# 一番左を足す
					if($hsp_end_hit[$i][$order[0]] < $subject_end[$i]){
						$hsp_stt[$i][$order[0]] = $hsp_stt[$i][$order[0]] - $subject_end[$i] + $hsp_end_hit[$i][$order[0]];
					}
					# 一番右を足す
					if($hsp_stt_hit[$i][$order[$orders-1]] > 1){
						$hsp_end[$i][$order[$orders-1]] = $hsp_end[$i][$order[$orders-1]] + $hsp_stt_hit[$i][$order[$orders-1]];
					}
					# ギャップを埋める
					for(my $j=0;$j<$orders-1;$j++){
						if($hsp_stt_hit[$i][$order[$j]] > $hsp_end_hit[$i][$order[$j+1]]){
							$hsp_end[$i][$order[$j]] = $hsp_end[$i][$order[$j]] + $hsp_stt_hit[$i][$order[$j]] - $hsp_end_hit[$i][$order[$j+1]];
						}
					}
				}
			}

			print STDERR 'Localblast: querylength= '.$query_length."\n";
			draw_align();

			# 散布図用のデータを揃える。ややこしいのでアラインメント用とは別に。FPKM値はRSEMデータが必要なので別にやる。検索も重いので別が良い。

			$hitcontigs = 0;
			$hitcontig_score_max = 0;
			$hitcontig_length_max = 0;
			$hitcontig_score_min = 5000;
			$hitcontig_length_min = 5000;
			my $searchio2 = new Bio::SearchIO (-format => 'blast', -file => $blast_file);
			my $blastresult2 = $searchio2 -> next_result;
			while(my $hit = $blastresult2 -> next_hit){
				my $name = $hit -> name;
				$name =~ s/\|/_/;
				$hitcontig_name[$hitcontigs] = $name;
				my $ev = $hit -> significance;
				if($ev == 0){$ev = 1e-181;}
				$ev = (-log($ev)/log(10));
				$hitcontig_evalue[$hitcontigs] = $ev;
				$hitcontig_score[$hitcontigs] = $hit -> score;
				$hitcontig_length[$hitcontigs] = $hit -> length;
				if($hitcontig_score_max < $hitcontig_score[$hitcontigs] ){$hitcontig_score_max = $hitcontig_score[$hitcontigs];}
				if($hitcontig_length_max < $hitcontig_length[$hitcontigs] ){$hitcontig_length_max = $hitcontig_length[$hitcontigs];}
				if($hitcontig_score_min > $hitcontig_score[$hitcontigs] ){$hitcontig_score_min = $hitcontig_score[$hitcontigs];}
				if($hitcontig_length_min > $hitcontig_length[$hitcontigs] ){$hitcontig_length_min = $hitcontig_length[$hitcontigs];}
				$hitcontigs++;
			}
			$hitcontig_resetflag = 1;
			blast2fasta();
			draw_scat();
		}else{
			print STDERR 'File is not a BLAST file: '.$buffer."\n";
		}
	}else{
		print STDERR 'File Not Found: '.$buffer."\n";
	}
}

#########
# NCBI BLASTファイルを選んでロードするボタンの動作
#########
sub findfile_ncbiblast(){
	my $buffer = $control_window -> getOpenFile( -filetypes => [['BLAST result','.txt'],['All files','*'] ]);
	if( -f $buffer){
		open IN, $buffer;my $test = <IN>;close IN;
		if($test =~ /BLAST/){
			$ncbiblast_file = $buffer;
			$ncbiblast_file_status = 'available';

			# データを変数に格納する
			@ncbi_chromat = ();
			$ncbi_chromat_max = 1;
			my $cutoff = 10;
			my $searchio_ncbi = new Bio::SearchIO (-format => 'blast', -file => $ncbiblast_file);
			my $blastresult_ncbi = $searchio_ncbi -> next_result;
			NCBI:while(my $hit = $blastresult_ncbi -> next_hit){
				my $judge = $hit -> significance;
				if($judge > 0){$judge = (-log($judge));}
				if($judge == 0){$judge = 181;}
				if($judge < $cutoff){last NCBI;}
				while (my $hsp = $hit -> next_hsp){
					my $start = $hsp -> start('query');
					my $end = $hsp -> end('query');
					my $eval = $hsp -> evalue;
					if($eval > 0){$eval = (-log($eval));}
					if($eval == 0){$eval = 181;}
					$ncbi_chromat[$start] += $eval;	if($ncbi_chromat_max < $ncbi_chromat[$start]){$ncbi_chromat_max = $ncbi_chromat[$start];}
					$ncbi_chromat[$end] += $eval;	if($ncbi_chromat_max < $ncbi_chromat[$end]  ){$ncbi_chromat_max = $ncbi_chromat[$end];  }
				}
			}
			print STDERR 'NCBI blast: maxpeak = '.$ncbi_chromat_max."\n";

		}else{
			print STDERR 'File is not a BLAST file: '.$buffer."\n";
		}
	}else{
		print STDERR 'File Not Found: '.$buffer."\n";
	}
}

#########
# RSEM 結果ファイルのパスをユーザーに指定させる
#########
sub findfile_rsem(){
	$rsem_file = $control_window -> getOpenFile( -filetypes => [['RSEM result','.txt'],['All files','*'] ]);
	# 一応チェック
	if( -f $rsem_file){
		open FILE, $rsem_file;
		my $line = <FILE>;
		if($line =~ /FPKM|TPM|expression/){$rsem_file_status = 'available';}else{$rsem_file_status = 'invalid';}
		close FILE;
		if($rsem_file_status eq 'available'){load_rsem();}
	}
}

#########
# Trinity 結果ファイルのパスをユーザーに指定させる
#########
sub findfile_trinity(){
	$trinity_file = $control_window -> getOpenFile( -filetypes => [['Trinity output','.fasta'],['All files','*']]);
	# 一応チェック
	if( -f $trinity_file){
		open FILE, $trinity_file;
		my $line = <FILE>;
		if($line =~ /\>/){$trinity_file_status = 'available';}else{$trinity_file_status = 'invalid';}
		close FILE;
		if($trinity_file_status eq 'available'){load_trinity();}
	}
}


#########
# 現在の設定をセーブ
#########
sub save_settings(){
	open OUT, '>RNAseqViewer_settings.txt';		# 前回のセッティング情報を覚えておくようにセット
		print OUT $rsem_file."\n";
		print OUT $trinity_file."\n";
		print OUT $local_library."\n";
		for(my $i=0;$i<$scats;$i++){
			print OUT $scat_xfpkm[$i]."\t".$scat_yfpkm[$i]."\t".$scat_label[$i]."\n";
		}
	close OUT;
}

#########
# RSEM結果とTrinityデータをロード。ついでにBLAST から得たname を元にlength をロード
#########
sub load_data(){
	load_rsem();
	load_trinity();
	save_settings();
}
sub load_rsem(){
	$fpkms = 0;
	$trinities = 0;
	$contigs = 0;
	$reset_scat_back = 1;

	open RSM, $rsem_file;
		# 最初の行が邪魔なのでスキップ。ついでにFPKM値が何サンプル分あるか調べる
		$firstlinersem = <RSM>;
		my @firstlinersem_split = split(/\t/,$firstlinersem);
		FST:for(my $i=2;$i<100;$i++){
			if($firstlinersem_split[$i] =~ /FPKM|TPM|expression/){
				$sample_name[$fpkms] = $firstlinersem_split[$i];
				$sample_name[$fpkms] =~ s/FPKM|TPM|expression//g;
				$fpkms++;
			}else{
				last FST;
			}
		}
		$scatx_listbox -> delete(0,'end');
		$scaty_listbox -> delete(0,'end');
		for(my $i=0;$i<$fpkms;$i++){
			$scatx_listbox -> insert($i,$i);
			$scaty_listbox -> insert($i,$i);
		}
		for(my $sc=0;$sc<$scats;$sc++){
			if($scat_xfpkm[$sc] >= $fpkms){$scat_xfpkm[$sc]=0;}
			if($scat_yfpkm[$sc] >= $fpkms){$scat_yfpkm[$sc]=0;}
		}
		$scat_listbox -> delete(0,'end');
		for(my $i=0;$i<$scats;$i++){$scat_listbox -> insert($i, 'x#'.$scat_xfpkm[$i].', y#'.$scat_yfpkm[$i].', '.$scat_label[$i]);}
		$reset_scat_back = 1;

		print STDERR 'RNAseq samples: '.$fpkms."\n";
		print STDERR 'Loading RSEM result...';

	while(my $line = <RSM>){	# ロード。計算の邪魔な|とかFPKM=0 とかをなんとかしておく
		chomp $line;
		my @s = split(/\t/,$line);
		$s[0] =~ s/\|/_/;							# 縦棒は計算式の邪魔になるので _ で置換する。
		$contig_line[$contigs] = $line;
		$contig_name[$contigs] = $s[0];
		$contig_length[$contigs] = $s[1];
		for(my $i=2;$i<2+$fpkms;$i++){
			if($s[$i] == 0){$s[$i] = 0.001;}		# FPKM = 0 だとlog計算できないので最低値のひとケタ下の値 0.001 を代わりに入れる
			$contig_fpkm[$contigs][$i-2] = $s[$i];
		}
		$contigs++;
		if($contigs % 10000 == 0){print STDERR '.';}
	}
	close RSM;
	print STDERR $contigs."\n";
	$rsem_file_status = $contigs. ' entries';
	draw_scat();
}
sub load_trinity(){
	$trinities = 0;

	# Trinity.fasta をロード。
	print STDERR 'Loading Trinity.fasta...';
	my $seqin = Bio::SeqIO->new(-file => $trinity_file, -format => 'fasta');
	while(my $seqobj = $seqin -> next_seq){
		$trinity_name[$trinities] = $seqobj -> id;
		$trinity_name[$trinities] =~ s/\|/_/;
		$trinity_seq[$trinities] = $seqobj -> seq;
		$trinities++;
		if($trinities % 10000 == 0){print STDERR '.';}
	}
	print STDERR $trinities."\n";
	$trinity_file_status = $trinities. ' entries';
}


#########
# ふたつのblast結果を元にアラインメント図を描く
#########

sub draw_align(){

	my $z = $align_zoom;	# GDのリサイズコマンドが何故か正常に動作しないので描画段階で拡大縮小することにした

	my $leftedge = 0;
	for(my $i=0;$i<$ranks;$i++){if($leftedge > $hsp_stt[$i][$ord[$i][0]]){$leftedge = $hsp_stt[$i][$ord[$i][0]];}}
	if($leftedge < (-400) ){$leftedge = (-400);}
	my $rightedge = $query_length;
	for(my $i=0;$i<$ranks;$i++){if($rightedge < $hsp_end[$i][$ord[$i][$ords[$i]-1]]){$rightedge = $hsp_end[$i][$ord[$i][$ords[$i]-1]];}}
	if($rightedge < $query_length+100 ){$rightedge = $query_length+100;}

	$leftedge = $leftedge * $z;
	$rightedge = $rightedge * $z;

	my $upmergin = 200 * $z;
	my $downmergin = 50 * $z;
	my $leftmergin = 400 * $z;
	my $rightmergin = 50 * $z;

	my $width = 60 * $z;
	my $bar_width = int($width * 0.65);
	my $fill_width = int($width * 0.5);
	my $overap_width = int($width * 0.15);
	my $chromat_baseline = int($width * 1.2);
	my $query_width = int($width * 1);

	my $graphsize_x=$rightedge - $leftedge + $leftmergin + $rightmergin;
	my $graphsize_y=$ranks * $width + $upmergin + $downmergin;
	my $image_align = GD::Image -> new($graphsize_x, $graphsize_y);
	my $white = $image_align -> colorAllocate(255,255,255);
	my $black = $image_align -> colorAllocate(0,0,0);
	my $red = $image_align -> colorAllocate(255,0,0);
	my $blue = $image_align -> colorAllocate(0,127,255);
	my $green = $image_align -> colorAllocate(0,127,0);

	if(Exists($alignment_window)){}else{
		$alignment_window = $control_window -> Toplevel();
		$alignment_window -> title('Alignment View '.$workidentifier);
		$alignment_window -> geometry("1000x700");
		$align_frame = $alignment_window -> Scrolled('Canvas',-scrollbars => 'se') -> pack(-fill => 'both', -expand => 1);
		$align_canvas = $align_frame -> Subwidget('scrolled') -> pack(-fill => 'both', -expand => 1);
	}
	$align_canvas -> delete('all');

	# NCBI blast でhspのエッジの位置をevalで量化したクロマトグラム
	if($ncbi_chromat_max > 1){
		for(my $i=0;$i<$query_length-1;$i++){
			my $x1 = $i*$z + $leftmergin - $leftedge;
			my $x2 = $i*$z + $leftmergin - $leftedge +1 ;
			my $y1 = $upmergin - $chromat_baseline - $ncbi_chromat[$i]/$ncbi_chromat_max*$upmergin;
			my $y2 = $upmergin - $chromat_baseline - $ncbi_chromat[$i+1]/$ncbi_chromat_max*$upmergin;
			$image_align -> line($x1,$y1,$x2,$y2,$red);
			$align_canvas -> create('line', $x1,$y1,$x2,$y2, -fill => 'red', -width => 2);
		}
		my $pointsize = 24 * $z;
		$image_align -> stringFT($red,$ttfont,$pointsize,     0,  50 * $z,  $upmergin - $width*2 + $fill_width, 'blastx edge');
		$align_canvas -> create('text', 50 * $z + 12*$pointsize/4,  $upmergin - $width*2 - $width/2 + $fill_width, -text => 'blastx edge',  -font => 'Arial '.$pointsize);
 	}
	# アラインメントの絵
	if($query_length > 1){
		for(my $i=0;$i<$ranks;$i++){
				my $prev_x2 = -1000;
				my $overap_switch = 0;
			for(my $j=0;$j<$ords[$i];$j++){
				if($hsp_stt[$i][$ord[$i][$j]] < $prev_x2 ){$overap_switch = $overap_width - $overap_switch;}
				my $x1 = $leftmergin - $leftedge + $hsp_stt[$i][$ord[$i][$j]] * $z;
				my $x2 = $leftmergin - $leftedge + $hsp_end[$i][$ord[$i][$j]] * $z;
				my $y1 = $upmergin + $overap_switch + $i*$width;
				my $y2 = $upmergin + $overap_switch + $i*$width + $bar_width;
				$prev_x2 = $hsp_end[$i][$ord[$i][$j]];
				$image_align -> filledRectangle($x1+1,$y1+1,$x2-1,$y2-1,$white);
				$align_canvas -> create('rectangle', $x1+1,$y1+1,$x2-1,$y2-1, -fill => 'white', -outline => 'black', -width => 2);
				for(my $k=0;$k<$cons[$i][$ord[$i][$j]];$k++){
					my $st = $leftmergin - $leftedge + $cons_stt[$i][$ord[$i][$j]][$k] * $z;
					my $en = $leftmergin - $leftedge + $cons_end[$i][$ord[$i][$j]][$k] * $z;
					my $hi = $upmergin + $overap_switch + $i*$width;
					my $lo = $upmergin + $overap_switch + $i*$width + $fill_width;
					$image_align ->       filledRectangle($st,$hi,$en,$lo, $blue);
					$align_canvas -> create('rectangle', $st,$hi,$en,$lo, -fill => 'blue', -outline => 'blue', -width => 1);
					if($j>0 && $k>0){
						$image_align ->       line(         $st-1,$hi+$fill_width*4/5,$st-1,$lo, $white);
						$align_canvas -> create('line', $st-1,$hi+$fill_width*4/5,$st-1,$lo, -fill => 'white', -width => 1);
					}	# 縮小するとギャップが潰れてしまうので左端の１ピクセル前に白い線を入れる
				}

				my $x1 = $leftmergin - $leftedge + $hsp_stt[$i][$ord[$i][$j]] * $z;
				my $x2 = $leftmergin - $leftedge + $hsp_end[$i][$ord[$i][$j]] * $z;
				my $y1 = $upmergin + $overap_switch + $i*$width;
				my $y2 = $upmergin + $overap_switch + $i*$width + $bar_width;
				$prev_x2 = $hsp_end[$i][$ord[$i][$j]];
				$image_align -> rectangle($x1,$y1,$x2,$y2,$black);
				$image_align -> rectangle($x1-1,$y1-1,$x2+1,$y2+1,$black);
				$align_canvas -> create('rectangle', $x1+1,$y1+1,$x2-1,$y2-1, -outline => 'black', -width => 2);
			}
		}

		# 文字を入れる
		my $pointsize = 24 * $z;
		for(my $i=0;$i<$ranks;$i++){
			my $label = $hit_name[$i];
			if ($subject_strand[$i] == (-1)){$label = $label.' (-)';}
			my $labelx =  50 * $z;
			my $labely =  $upmergin + $i*$width + $bar_width;
			$image_align -> stringFT($black,$ttfont,$pointsize,     0,$labelx,$labely,$label);
			$align_canvas -> create('text', $labelx+length($label)*$pointsize/4,$labely-$bar_width/2, -text => $label,  -font => 'Arial '.$pointsize);
		}

		# query のバーと文字を入れる
		my $x1 = $leftmergin - $leftedge;
		my $x2 = $leftmergin - $leftedge + $query_length * $z;
		my $y1 = $upmergin  - $query_width;
		my $y2 = $upmergin  - $query_width;
		$image_align -> stringFT($green,$ttfont,$pointsize,     0,  50 * $z,  $upmergin - $width + $fill_width, 'Query');
		$image_align -> stringFT($black,$ttfont,$pointsize,     0,  $x1 - $pointsize ,  $y1 + $fill_width, '1'  );
		$image_align -> stringFT($black,$ttfont,$pointsize,     0,  $x2 + $pointsize/4, $y2 + $fill_width, $query_length);
		$align_canvas -> create('text', 50 * $z + 5*$pointsize/4 ,  $upmergin - $width + $fill_width - $bar_width/2, -text => 'Query',      -font => 'Arial '.$pointsize);
		$align_canvas -> create('text',  $x1 - $pointsize + $pointsize/4 ,    $y1 + $fill_width - $bar_width/2, -text => '1',                    -font => 'Arial '.$pointsize  );
		$align_canvas -> create('text',  $x2 + (length($query_length)+2)*$pointsize/4, $y2 + $fill_width - $bar_width/2, -text => $query_length, -font => 'Arial '.$pointsize);
		$image_align -> filledRectangle($x1,$y1,$x2,$y2+ $fill_width,$green);
		$image_align -> rectangle($x1,$y1,$x2,$y2+ $bar_width,$black);
		$image_align -> rectangle($x1-1,$y1-1,$x2+1,$y2+ $bar_width+1,$black);
		$align_canvas -> create('rectangle', $x1,$y1,$x2,$y2+ $fill_width, -fill => 'green');
		$align_canvas -> create('rectangle',$x1,$y1,$x2,$y2+ $bar_width,- outline => 'black', -width => 2);

#		my $png = $alignment_window -> Photo(-data => encode_base64($image_align -> png), -format => 'png');	# macでエラーを吐くので中止
		$align_canvas -> configure (-width => $graphsize_x, -height => $graphsize_y);	# これだけでは何故かサイズを小さい方には変更できない。
		$align_canvas -> configure (-scrollregion => [0,0,$graphsize_x,$graphsize_y]);
#		$align_canvas -> create('image', $graphsize_x/2, $graphsize_y/2, -image => $png);
#		$align_canvas -> create('line', 0, $graphsize_y, $graphsize_x, 0, -fill => $blue, -width => 2); # GDの像を表示できないので代わりに線を引く
		$alignment_window -> maxsize($graphsize_x+20,$graphsize_y+20);	# これで一応小さい方にも変更できるが、一回ウィンドウに触ってもらう必要がある。仕方ないからこれでいい
		my $align_savebutton_rec = $align_canvas -> create('rectangle', 4,4,4+48,4+24, -outline => 'black', -fill => 'gray');
		   $align_canvas -> bind($align_savebutton_rec, "<Button-1>" => [\&save_align]);
		my $align_savebutton_txt = $align_canvas -> create('text', (4+48)/2,(4+24)/2, -anchor => 'c', -text => 'Save');
		   $align_canvas -> bind($align_savebutton_txt, "<Button-1>" => [\&save_align]);

#		open OUT, ">temp_drawblast.png";
#		binmode OUT;
#		print OUT $image_align -> png();
#		close OUT;

		if(defined($align_savefile)){
			if($align_savefile =~ /\.png/){}else{$align_savefile = $align_savefile.'.png';}
			open OUT, '>'.$align_savefile;
			binmode OUT;
			print OUT $image_align -> png();
			close OUT;
			print STDERR 'save done: '.$align_savefile."\n";
			undef($align_savefile);
		}
	}
}

#########
# アラインメント図をセーブする。
#########

sub save_align(){
	$align_savefile = $alignment_window -> getSaveFile( -filetypes => [['PNG image','.png'],['All files','*']]);
	if(defined($align_savefile)){
		draw_align();
	}else{
		$align_savefile ='RNAseqViewer_blastalign.png';
		print STDERR 'save canceled'."\n";
	}
}

#########
# 散布図の背景を描く。毎回描くと遅いので必要な時だけ呼ぶ
#########

sub draw_scat_background(){
	if($reset_scat_back == 1){if($contigs>0){if($scats>0){

		my $z = $scat_zoom;
		$graphsize=900*$z;	# キャンバスいっこ分のサイズ
		$backimage = GD::Image ->new($graphsize * $scats,$graphsize);	# 設定数分の枠を用意する
		my $white = $backimage -> colorAllocate(255,255,255);
		my $black = $backimage -> colorAllocate(0,0,0);
		my $red = $backimage -> colorAllocate(255,0,0);
		my $green = $backimage -> colorAllocate(0,127,0);
		my $blue = $backimage -> colorAllocate(127,127,255);
		my $darkblue = $backimage -> colorAllocate(0,0,255);
		my $gray = $backimage -> colorAllocate(127,127,127);
		#$backimage -> transparent($white);
		$backimage -> interlaced('true');
		my $mergin = 80*$z;		# 左と下の余白
		my $scale = 100*$z;		# 10倍目盛のサイズ

		for(my $sc =0; $sc<$scats; $sc++){

			# 目盛線と枠を表示
			my $xzero = $mergin + $graphsize * $sc;
			my $yzero = $graphsize-$mergin;
			{
				my $fontsize = 16*$z;
				for(my $i=1;$i<8;$i++){
					$backimage -> line($xzero+$i*$scale,     $yzero,          $xzero+$i*$scale,      $yzero-8*$scale, $gray);		# 横軸補助目盛線
					$backimage -> line($xzero+$i*$scale,     $yzero,          $xzero+$i*$scale,      $yzero-8,             $black);		# 横軸めもり
					$backimage -> line($xzero+1+$i*$scale, $yzero,          $xzero+$i*$scale+1,  $yzero-8,             $black);		# 横軸めもり
					$backimage -> line($xzero,          $yzero-$i*$scale,     $xzero+8*$scale, $yzero-$i*$scale,      $gray);		# 縦軸補助目盛線
					$backimage -> line($xzero,          $yzero-$i*$scale,     $xzero+8            , $yzero-$i*$scale,      $black);		# 縦軸めもり
					$backimage -> line($xzero,          $yzero-$i*$scale-1,  $xzero+8            , $yzero-$i*$scale-1,  $black);		# 縦軸めもり
					my $label = 10**($i-3);
					if($scat_label[$sc] =~ /olcano/){$label=($i-4);}
					$backimage -> stringFT($black, $ttfont, $fontsize, 0, $xzero+4*$z+$i*$scale-length($label)*$fontsize/2,$yzero+8*$z+$fontsize,$label);		# 横軸の数字
					if($scat_label[$sc] =~ /olcano/){$label=$i;}
					$backimage -> stringFT($black, $ttfont, $fontsize, 0, $xzero-10*$z-length($label)*$fontsize/(1.5),$yzero-$i*$scale+$fontsize/2,$label);		# 縦軸の数字
				}
				$backimage -> line($xzero, $yzero, $xzero+8*$scale, $yzero,         $black);	# 横の目盛線
				$backimage -> line($xzero, $yzero, $xzero,          $yzero-8*$scale,$black);	# 縦の目盛線
			}
			# 軸ラベルを表示
			{	my $fontsize = 24*$z;
#				my $label = 'Sample #'.$scat_xfpkm[$sc].' (FPKM)';	# 横軸のラベル
				my $label = $sample_name[$scat_xfpkm[$sc]];			# 横軸のラベル
				   if($scat_label[$sc] =~ /olcano/){$label=$label.'(log4)';}
				$backimage -> stringFT($black, $ttfont, $fontsize, 0,      $xzero+8*$scale/2-length($label)*$fontsize/3,$yzero+$fontsize+40*$z,$label);
#				   $label = 'Sample #'.$scat_yfpkm[$sc].' (FPKM)';	# 縦軸のラベル
				   $label = $sample_name[$scat_yfpkm[$sc]];			# 縦軸のラベル
				   if($scat_label[$sc] =~ /olcano/){$label=$label.'(-log10)';}
				$backimage -> stringFT($black, $ttfont, $fontsize, 1.5708, $xzero-$fontsize-32*$z, $yzero-8*$scale/2+length($label)*$fontsize/3,$label);
				   $label = $scat_label[$sc];						# 図のラベル
				   if($scat_label[$sc] =~ /olcano/){$label='';}
				$backimage -> stringFT($black, $ttfont, $fontsize, 0,      $xzero+8*$scale-length($label)*$fontsize/1.5,$yzero-$fontsize,$label);
			}

			# 背景をプロットする
			my $pointsize = 2;
			for(my $i=0;$i<$contigs;$i++){
				my $x1 = $xzero + $scale * (log($contig_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
				my $y1 = $yzero - $scale * (log($contig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
				my $x2 = $xzero + $scale * (log($contig_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
				my $y2 = $yzero - $scale * (log($contig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
				if($scat_label[$sc] =~ /olcano/){
					$x1 = $xzero + $scale * (log($contig_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
					$y1 = $yzero - $scale * (-log($contig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
					$x2 = $xzero + $scale * (log($contig_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
					$y2 = $yzero - $scale * (-log($contig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
				}
				$backimage -> filledRectangle($x1,$y1,$x2,$y2,$blue);
			}
			
			# x=y の線を引く
			if($scat_label[$sc] =~ /olcano/){
				$backimage -> line($xzero+4*$scale,   $yzero,   $xzero+4*$scale,   $yzero-8*$scale,   $green);
				$backimage -> line($xzero+4*$scale+1, $yzero,   $xzero+4*$scale+1, $yzero-8*$scale,   $green);	# ちょっと細いので太くする
				$backimage -> line($xzero+4*$scale-1, $yzero,   $xzero+4*$scale-1, $yzero-8*$scale,   $green);
			}else{
				$backimage -> line($xzero,   $yzero,   $xzero+8*$scale,   $yzero-8*$scale,   $green);
				$backimage -> line($xzero+1, $yzero,   $xzero+8*$scale,   $yzero-8*$scale+1, $green);	# ちょっと細いので太くする
				$backimage -> line($xzero,   $yzero-1, $xzero+8*$scale-1, $yzero-8*$scale,   $green);
			}
		}
		$reset_scat_back = 0;
	}}}
}


#########
# 散布図を描く。背景は別に描いて高速化する
#########

sub draw_scat(){
	print STDERR 'call draw_scat'."\n";

			my $z = $scat_zoom;
			$graphsize=900*$z;	# キャンバスいっこ分のサイズ
			my $mergin = 80*$z;	# 左と下の余白
			my $scale = 100*$z;	# 10倍目盛のサイズ
			my $pointsize = 5;
			my @fieldx;
			my @fieldy;
			my $fields=0;

	if($contigs>0){if($scats>0){
		if($reset_scat_back == 1){draw_scat_background();}			# 背景画像のリセットフラグの参照が二重だが一応とっとく
		if($hitcontig_resetflag == 1){search_hitcontig_fpkm();}		# Local blast が更新された場合は描画前にfpkm値の検索をやり直す
		if($textwidget_refreshflag == 1){input_fasta();}			# スキャッタープロットの前にエントリーリストを確認する。
		if(Exists($plot_window)){}else{
			$plot_window = $control_window -> Toplevel();
			$plot_window -> title('ScatterPlot View '.$workidentifier);
			$plot_window -> geometry("1300x650");
			$plot_frame = $plot_window -> Scrolled('Canvas',-scrollbars => 'se') -> pack(-fill => 'both', -expand => 1);
			$plot_canvas = $plot_frame -> Subwidget('canvas') -> pack(-fill => 'both', -expand => 1);
		}
#		my $png = $plot_window -> Photo(-data => encode_base64($backimage -> png), -format => 'png'); # Photo が機能しないので
		$plot_canvas -> delete('all');
		$plot_frame -> configure (-width => $graphsize * $scats, -height => $graphsize);
		$plot_canvas -> configure (-width => $graphsize * $scats, -height => $graphsize);	# これだけでは何故かサイズを小さい方には変更できない。
		$plot_canvas -> configure (-scrollregion => [0,0,$graphsize * $scats,$graphsize]);
#		my $background_plotimage = $plot_canvas -> create('image', $graphsize * $scats/2, $graphsize/2, -image => $png);
		# GDの像を表示できないので代わりに枠線と斜め線を引く
		for(my $sc=0;$sc<$scats;$sc++){
			my $xzero = $mergin + $graphsize * $sc;
			my $yzero = $graphsize-$mergin;
			$plot_canvas -> create('rectangle', $xzero, $yzero, $xzero+8*$scale, $yzero-8*$scale, -fill => 'white', -outline => 'black',-width => 1);
				my $fontsize = 16*$z;
				for(my $i=1;$i<8;$i++){
					$plot_canvas -> create('line', $xzero+$i*$scale,$yzero,          $xzero+$i*$scale,$yzero-8*$scale, -fill => 'gray', -width => 1);		# 横軸補助目盛線
					$plot_canvas -> create('line', $xzero+$i*$scale,$yzero,          $xzero+$i*$scale,$yzero-8,        -fill => 'black', -width => 2);		# 横軸めもり
					$plot_canvas -> create('line', $xzero,          $yzero-$i*$scale,$xzero+8*$scale, $yzero-$i*$scale, -fill => 'gray', -width => 1);		# 縦軸補助目盛線
					$plot_canvas -> create('line', $xzero,          $yzero-$i*$scale,$xzero+8       , $yzero-$i*$scale, -fill => 'black', -width => 2);		# 縦軸めもり
					my $label = 10**($i-3);
					$plot_canvas -> create('text', $xzero+$i*$scale-length($label)*$fontsize/16,$yzero+$fontsize, -text => $label, -font => 'Arial '.$fontsize); # 横軸の数字
					$plot_canvas -> create('text', $xzero-14*$z-length($label)*$fontsize/4,$yzero-$i*$scale, -text => $label, -font => 'Arial '.$fontsize); # 縦軸の数字
				}
			if($scat_label[$sc] =~ /olcano/){
				$plot_canvas -> create('line', $xzero+4*$scale, $yzero, $xzero+4*$scale, $yzero-8*$scale, -fill => 'blue', -width => 2);
			}else{
				$plot_canvas -> create('line', $xzero, $yzero, $xzero+8*$scale, $yzero-8*$scale, -fill => 'blue', -width => 2);
			}
			$fontsize = 24*$z;
			my $label = 'Sample #'.$scat_yfpkm[$sc].' / #'.$scat_xfpkm[$sc].': '.$scat_label[$sc];						# 図のラベル
			$plot_canvas -> create('text', $xzero+7*$scale-length($label)*$fontsize/4,$yzero-$fontsize-$scale*0.3, -text => $label, -font => 'Arial '.$fontsize);
#			$backimage -> stringFT($black, $ttfont, $fontsize, 0,      $xzero+8*$scale-length($label)*$fontsize/1.5,$yzero-$fontsize,$label);
		}
		if($query_length == 0 && $keyword_hits == 0){$plot_canvas -> bind($background_plotimage, "<Button-1>" => [\&create_by_userclick, Ev('x'), Ev('y')]);}
		$plot_window -> maxsize($graphsize * $scats+20,$graphsize+20);	# これで一応小さい方にも変更できるが、一回ウィンドウに触ってもらう必要がある。仕方ないからこれでいい

		$saveimage_scat = $backimage -> clone();
		my $white = $saveimage_scat -> colorAllocate(255,255,255);
		my $black = $saveimage_scat -> colorAllocate(0,0,0);
		my $red = $saveimage_scat -> colorAllocate(255,0,0);
		my $savecolor_index;
		my @savecolor;
		for(my $i=0;$i<128;$i++){
			my $r =     int(511* $i /127);if($r>255){$r=255;}
			my $g = 511-int(511* $i /127);if($g>255){$g=255;}
			my $b = 0;
			$savecolor[$i] = $saveimage_scat -> colorAllocate($r,$g,$b);
		}


		# Local Blast の結果をプロットに反映させる

		if($query_length > 1){
			for(my $sc =0; $sc<$scats; $sc++){
				my $xzero = $mergin + $graphsize * $sc;
				my $yzero = $graphsize-$mergin;
				# プロットする
				for(my $i=$hitcontigs-1;$i>=0;$i--){
					my $x1 = $xzero + $scale * (log($hitcontig_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
					my $y1 = $yzero - $scale * (log($hitcontig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
					my $x2 = $xzero + $scale * (log($hitcontig_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
					my $y2 = $yzero - $scale * (log($hitcontig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
					if($scat_label[$sc] =~ /olcano/){
						$x1 = $xzero + $scale * (log($hitcontig_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
						$y1 = $yzero - $scale * (-log($hitcontig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
						$x2 = $xzero + $scale * (log($hitcontig_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
						$y2 = $yzero - $scale * (-log($hitcontig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
					}

					my $pointcolor;
					if($scat_plotcolor_switch == 0){
						my $r = int(511 * ($hitcontig_score[$i] - $hitcontig_score_min)/ ($hitcontig_score_max - $hitcontig_score_min));
						if($r>255){$r = 255;};$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g = 511-int(511 * ($hitcontig_score[$i] - $hitcontig_score_min) / ($hitcontig_score_max - $hitcontig_score_min));
						if($g>255){$g = 255;};$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b = sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = int(127*($hitcontig_score[$i] - $hitcontig_score_min)/ ($hitcontig_score_max - $hitcontig_score_min));
					}
					if($scat_plotcolor_switch == 1){
						my $r = int(511 * $hitcontig_evalue[$i] / 181);
						if($r>255){$r = 255;};$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g = 511-int(511 * $hitcontig_evalue[$i] / 181);
						if($g>255){$g = 255;};$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b = sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = int(127 * $hitcontig_evalue[$i] / 181);
					}
					if($scat_plotcolor_switch == 2){
						my $r = int(511 * ($hitcontig_length[$i] - $hitcontig_length_min)/ ($hitcontig_length_max - $hitcontig_length_min));
						if($r>255){$r = 255;};$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g = 511-int(511 * ($hitcontig_length[$i] - $hitcontig_length_min) / ($hitcontig_length_max - $hitcontig_length_min));
						if($g>255){$g = 255;};$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b = sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = int(127 * ($hitcontig_length[$i] - $hitcontig_length_min)/ ($hitcontig_length_max - $hitcontig_length_min));
					}
					my $plot = $plot_canvas -> create('rectangle', $x1,  $y1,  $x2,  $y2,  -outline => 'black', -fill => $pointcolor, -width => 2);
					$saveimage_scat -> filledRectangle($x1,$y1,$x2,$y2,$black);
					$saveimage_scat -> filledRectangle($x1+1,$y1+1,$x2-1,$y2-1,$savecolor[$savecolor_index]);
					$plot_canvas -> bind($plot, "<Button-1>" => [\&select_by_userclick, Ev('x'), Ev('y')]);
					$hitcontig_plotx[$sc][$i] = ($x1+$x2)/2;	# find が機能しないので仕方なく
					$hitcontig_ploty[$sc][$i] = ($y1+$y2)/2;

					# プロットした座標をとっとく
					$fieldx[$fields]=($x1+$x2)/2;
					$fieldy[$fields]=($y1+$y2)/2;
					$fields++;
				}
			}
		}

		if($query_length == 0){if($keyword_hits>0){		# キーワード検索のときの色
			for(my $sc =0; $sc<$scats; $sc++){
				my $xzero = $mergin + $graphsize * $sc;
				my $yzero = $graphsize-$mergin;
				# プロットする
				for(my $i=$keyword_hits-1;$i>=0;$i--){
					my $x1 = $xzero + $scale * (log($keyword_hit_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
					my $y1 = $yzero - $scale * (log($keyword_hit_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
					my $x2 = $xzero + $scale * (log($keyword_hit_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
					my $y2 = $yzero - $scale * (log($keyword_hit_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
					if($scat_label[$sc] =~ /olcano/){
						$x1 = $xzero + $scale * (log($keyword_hit_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
						$y1 = $yzero - $scale * (-log($keyword_hit_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
						$x2 = $xzero + $scale * (log($keyword_hit_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
						$y2 = $yzero - $scale * (-log($keyword_hit_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
					}

					my $pointcolor;
					if($keyword_hit_length_max - $keyword_hit_length_min > 0){		# ゼロ割エラー出ちゃうので。
						my $r = int(511 * ($keyword_hit_length[$i] - $keyword_hit_length_min)/ ($keyword_hit_length_max - $keyword_hit_length_min));
						if($r>255){$r = 255;};$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g = 511-int(511 * ($keyword_hit_length[$i] - $keyword_hit_length_min) / ($keyword_hit_length_max - $keyword_hit_length_min));
						if($g>255){$g = 255;};$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b = sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = int(127 * ($keyword_hit_length[$i] - $keyword_hit_length_min)/ ($keyword_hit_length_max - $keyword_hit_length_min));
#						$savecolor_index = int(127 * ($hitcontig_length[$i] - $hitcontig_length_min)/ ($hitcontig_length_max - $hitcontig_length_min));
					}else{
						my $r=255;$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g=255;$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b=sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = 63;
					}

					if($scat_plotcolor_switch == 0){	# キーワード検索のハイライトの場合に黄色一色にするスイッチとして使用してみる
						my $r=255;$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g=255;$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b=sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = 63;
					}

					my $plot = $plot_canvas -> create('rectangle', $x1,  $y1,  $x2,  $y2,  -outline => 'black', -fill => $pointcolor, -width => 2);
					$saveimage_scat -> filledRectangle($x1,$y1,$x2,$y2,$black);
					$saveimage_scat -> filledRectangle($x1+2,$y1+2,$x2-2,$y2-2,$savecolor[$savecolor_index]);
					$plot_canvas -> bind($plot, "<Button-1>" => [\&select_by_userclick, Ev('x'), Ev('y')]);
					$keyword_hit_plotx[$sc][$i] = ($x1+$x2)/2;	# find が機能しないので仕方なく
					$keyword_hit_ploty[$sc][$i] = ($y1+$y2)/2;

					# プロットした座標をとっとく
					$fieldx[$fields]=($x1+$x2)/2;
					$fieldy[$fields]=($y1+$y2)/2;
					$fields++;
				}
				
				# キーワードを表示する
				my $label = $keyword_input.' ('.$keyword_hits.' genes)';
				$label =~ s/\|/\,\ /g;			# or 演算子は検索に有効だけど見づらいのでコンマと空白を入れる
				my $title_x = $xzero+$scale*1;	# 左から１マス目
				my $title_y = $yzero-$scale*6;	# 下から６マス目
				my $fontsize = 24*$z;			# 軸ラベルと同じサイズ= 24 * $z
				$saveimage_scat -> stringFT($black, $ttfont, $fontsize, 0, $title_x, $title_y, $label);
			}
		}}

		if($fastas > 0){
			for(my $sc =0; $sc<$scats; $sc++){
				my $xzero = $mergin + $graphsize * $sc;
				my $yzero = $graphsize-$mergin;
				# ユーザーの指定コンティグをハイライト表示にする
					for(my $i=$fastas-1;$i>=0;$i--){
						if(defined($fasta_fpkm[$i][$scat_xfpkm[$sc]])){
							my $x1 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
							my $y1 = $yzero - $scale * (log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
							my $x2 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
							my $y2 = $yzero - $scale * (log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
							if($scat_label[$sc] =~ /olcano/){
								$x1 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
								$y1 = $yzero - $scale * (-log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
								$x2 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
								$y2 = $yzero - $scale * (-log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
							}
							$plot_canvas -> create('rectangle', $x1,  $y1,  $x2,  $y2,  -outline => 'black', -width => 2);
							$saveimage_scat -> rectangle($x1,$y1,$x2,$y2,$black);
							$saveimage_scat -> rectangle($x1+2,$y1+2,$x2-2,$y2-2,$black);
							$x1 = $x1 - $pointsize;$y1 = $y1 - $pointsize;
							$x2 = $x2 + $pointsize;$y2 = $y2 + $pointsize;
							$plot_canvas -> create('oval',      $x1,  $y1,  $x2,  $y2,  -outline => 'red', -width => 2);
							$saveimage_scat -> arc(($x1+$x2)/2, ($y1+$y2)/2, $pointsize*4, $pointsize*4, 0, 360, $red);
							$fieldx[$fields]=($x1+$x2)/2;
							$fieldy[$fields]=($y1+$y2)/2;
							$fields++;
						}
					}
				my $centerx=0;
				my $centery=0;
				my $centerc=0;
				if($fields>0){
					# 文字をずらす方向を決めるために今あるスポット中心あたりの座標を探す
					for(my $i=0;$i<$fields;$i++){
						if($fieldx[$i]>$xzero){if($fieldy[$i]<$yzero){
							$centerx = $centerx + $fieldx[$i];
							$centery = $centery + $fieldy[$i];
							$centerc++;
						}}
					}
					if($centerc>0){
						$centerx = $centerx / $centerc;
						$centery = $centery / $centerc;
					}else{
						$centerx = $xzero + $scale * 4;
						$centery = $yzero - $scale * 4;
					}
				}else{
					$centerx = $xzero + $scale * 4;
					$centery = $yzero - $scale * 4;
				}
					for(my $i=$fastas-1;$i>=0;$i--){
						if(defined($fasta_fpkm[$i][$scat_xfpkm[$sc]])){
							if($scat_hitextswitch == 1){
							if($i<$scat_showtext_limit){

								# 文字が目立つ位置まで移動する
								my $x1 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
								my $y1 = $yzero - $scale * (log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
								my $x2 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
								my $y2 = $yzero - $scale * (log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
								if($scat_label[$sc] =~ /olcano/){
									$x1 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
									$y1 = $yzero - $scale * (-log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
									$x2 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
									$y2 = $yzero - $scale * (-log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
								}
								my $posx = ($x1+$x2)/2;
								my $posy = ($y1+$y2)/2;
								my $effort =0;
								my $nearest=0;
								DST:while($nearest<150*$z){
									# 良い感じのコンティグの場合
									if($fasta_fpkm[$i][$scat_xfpkm[$sc]] >0.001 && $fasta_fpkm[$i][$scat_yfpkm[$sc]] >0.001){
										my $vector = sqrt(($posx - $centerx)**2 + ($posy - $centery)**2);	# スポット群の中心点から遠ざかる方向へのベクトルの長さ
										my $force = 0;
										if($fasta_fpkm[$i][$scat_xfpkm[$sc]] - $fasta_fpkm[$i][$scat_yfpkm[$sc]] != 0){			# 右下か左上へのベクトル
											$force = ($fasta_fpkm[$i][$scat_xfpkm[$sc]] - $fasta_fpkm[$i][$scat_yfpkm[$sc]])
											    / abs($fasta_fpkm[$i][$scat_xfpkm[$sc]] - $fasta_fpkm[$i][$scat_yfpkm[$sc]]);
										}
										$posx = $posx + (10 * ($posx - $centerx)/$vector + rand(40)-20 + 6*$force)*$z;		# x方向のランダムを気持ち多め
										$posy = $posy + (10 * ($posy - $centery)/$vector + rand(30)-15 + 6*$force)*$z;
										if($posx < $xzero + 80*$z){$posx = ($x1+$x2)/2;$posy = ($y1+$y2)/2;}				# ハミ出たら戻す
										if($posy > $yzero - 12*$z){$posy = ($y1+$y2)/2;$posx = ($x1+$x2)/2;}
										if($posx > $xzero + $scale * 8 - 80*$z){$posx = $xzero + $scale * 8 - 80*$z;}		# ハミ出たら戻す
										if($posy < $yzero - $scale * 8 + 12*$z){$posy = $yzero - $scale * 8 + 12*$z;}
									# どっちかの軸がゼロになっている変なコンティグの場合
									}else{
										if($fasta_fpkm[$i][$scat_xfpkm[$sc]] == 0.001 && $fasta_fpkm[$i][$scat_yfpkm[$sc]] == 0.001){
											$posx = $posx + rand(10)*$z;$posy = $posy - rand(10)*$z;	# 両軸ともゼロの場合は右上方向に移動
										}else{
											if($fasta_fpkm[$i][$scat_xfpkm[$sc]] == 0.001){
												$posx = $posx + rand(10)*$z;						# x が初めからゼロの場合は右方向に移動
												$posy = $posy + (rand(20)-10)*$z;
												if($posx > $xzero + $scale * 2 - 80*$z){$posx = ($x1+$x2)/2;}		# ハミ出たら戻す
											}
											if($fasta_fpkm[$i][$scat_yfpkm[$sc]] == 0.001){
												$posx = $posx + (rand(20)-10)*$z;
												$posy = $posy - rand(10)*$z;
												if($posy < $yzero - $scale * 1 + 12*$z){$posy = ($y1+$y2)/2;}
											}
										}

									}

									# 一番近いプロットまでの距離を雑に求めて評価用の変数に入れる
									my $dist_min=1000;
									for(my $j=0;$j<$fields;$j++){
										my $dist = 80*$z+sqrt(abs($posx - $fieldx[$j])**2 + abs($posy - $fieldy[$j])**2);
										if($dist_min > $dist){$dist_min = $dist;}
									}
									$nearest = $dist_min;
									$effort++;
									if($effort>1000){last DST;}
								}
								my $text = $fasta_name[$i]; $text =~ s/_//g;
								$plot_canvas -> create('line', $posx,$posy, ($x1+$x2)/2, ($y1+$y2)/2, -width => 2, -fill => 'red');
					#			$plot_canvas -> create('rectangle', $posx-80*$z,$posy-12*$z,$posx+80*$z,$posy+12*$z, -width => 2, -fill => 'white', -outline => 'red');
								$plot_canvas -> create('rectangle', $posx-80,$posy-12,$posx+80,$posy+12, -width => 2, -fill => 'white', -outline => 'red');
								$saveimage_scat -> line($posx,$posy, ($x1+$x2)/2, ($y1+$y2)/2,$red);
					#			$saveimage_scat -> filledRectangle($posx-80*$z-1,$posy-12*$z-1,$posx+80*$z+1,$posy+12*$z+1,$red);
					#			$saveimage_scat -> filledRectangle($posx-80*$z+1,$posy-12*$z+1,$posx+80*$z-1,$posy+12*$z-1,$white);
								$saveimage_scat -> filledRectangle($posx-80-1,$posy-12-1,$posx+80+1,$posy+12+1,$red);
								$saveimage_scat -> filledRectangle($posx-80+1,$posy-12+1,$posx+80-1,$posy+12-1,$white);
					#			my $fontsize = 14*$z;
								my $fontsize = 14;
								$plot_canvas -> create('text',$posx,$posy, -text => $text, -font => 'Arial '.$fontsize);
								$saveimage_scat -> stringFT($black, $ttfont, $fontsize, 0, $posx-length($text)*$fontsize/3,$posy+$fontsize/2, $text);
								$fieldx[$fields]=$posx;
								$fieldy[$fields]=$posy;
								$fields++;
							}}
						}
					}
				#}
			}
		}

		# セーブボタンを表示
		my $plot_savebutton_rec = $plot_canvas -> create('rectangle', 4,4,4+48,4+24, -outline => 'black', -fill => 'gray');
		   $plot_canvas -> bind($plot_savebutton_rec, "<Button-1>" => [\&save_scat]);
		my $plot_savebutton_txt = $plot_canvas -> create('text', (4+48)/2,(4+24)/2, -anchor => 'c', -text => 'Save');
		   $plot_canvas -> bind($plot_savebutton_txt, "<Button-1>" => [\&save_scat]);

#		open OUT, ">temp_drawplot.png";
#		binmode OUT;
#		print OUT $saveimage_scat -> png();
#		close OUT;
#
		if(defined($scat_savefile)){
			if($scat_savefile =~ /\.png/){}else{$scat_savefile = $scat_savefile.'.png';}
			open OUT, '>'.$scat_savefile;
			binmode OUT;
			print OUT $saveimage_scat -> png();
			close OUT;
			print STDERR 'save done: '.$scat_savefile."\n";
#			undef($scat_savefile);
		}

	}}
}

#########
# Local blast でヒットしたコンティグのfpkm値を検索する
#########

sub search_hitcontig_fpkm(){
	my $count=0;
	for(my $j=0;$j<$hitcontigs;$j++){
		HIT:for(my $i=0;$i<$contigs;$i++){
			if($hitcontig_name[$j] eq $contig_name[$i]){
				for(my $k=0;$k<$fpkms;$k++){
					$hitcontig_fpkm[$j][$k] = $contig_fpkm[$i][$k];
				}
				$count++;
				last HIT;
			}
		}
	}
	print STDERR 'search: '.$count.' / '.$hitcontigs."\n";
	$hitcontig_resetflag = 0;
}

#########
# スキャッタープロット上でハイライトプロットがクリックされた時の動作。
#########

sub select_by_userclick(){
	my @argument = @_;	# 0=plotcanvasのobject, 1=x, 2=y が渡されているはず
	my $x = $argument[1];
	my $y = $argument[2];
	$x = $plot_canvas -> canvasx($x);
	$y = $plot_canvas -> canvasy($y);

	# 全プロットの中から一番近いやつを探す。できれば一番手前に表示されてるやつがいい
	my $near=10000;
	my $id;

	if($textwidget_refreshflag == 1){input_fasta();}

	if($query_length >0){
		for(my $sc=0;$sc<$scats;$sc++){
			SC:for(my $i=0;$i<$hitcontigs;$i++){
				my $dist = sqrt(($hitcontig_plotx[$sc][$i]-$x)**2 + ($hitcontig_ploty[$sc][$i]-$y)**2);
				if($near > $dist){
						$near = $dist;
						$id =$i;
				}
			}
		}
		print STDERR 'local blast hit#= '.$hitcontig_name[$id]."\n";
#		for(my $i=$fastas;$i>0;$i--){
#			$fasta_name[$i] = $fasta_name[$i-1];
#			$fasta_seq[$i] = $fasta_seq[$i-1];
#		}
#		$fasta_name[0] = $hitcontig_name[$id];
#		$fasta_seq[0] = $hitcontig_seq[$id];
		$fastas++;
		unshift(@fasta_name, $hitcontig_name[$id]);
		unshift(@fasta_seq,  $hitcontig_seq[$id]);
	}
	if($query_length == 0){if($keyword_hits>0){
		for(my $sc=0;$sc<$scats;$sc++){
			SC:for(my $i=0;$i<$keyword_hits;$i++){
				my $dist = sqrt(($keyword_hit_plotx[$sc][$i]-$x)**2 + ($keyword_hit_ploty[$sc][$i]-$y)**2);
				if($near > $dist){
						$near = $dist;
						$id =$i;
				}
			}
		}
		print STDERR 'keyword hit selected: # '.$keyword_hit_name[$id]."\n";
#		for(my $i=$fastas;$i>0;$i--){
#			$fasta_name[$i] = $fasta_name[$i-1];
#			$fasta_seq[$i] = $fasta_seq[$i-1];
#		}
#		$fasta_name[0] = $keyword_hit_name[$id];
#		$fasta_seq[0] = $keyword_hit_seq[$id];
		$fastas++;
		unshift(@fasta_name, $keyword_hit_name[$id]);
		unshift(@fasta_seq,  $keyword_hit_seq[$id]);
	}}

#	for(my $i=$fastas;$i>0;$i--)
#		for(my $j=0;$j<$fpkms;$j++){
#			$fasta_fpkm[$i][$j]=$fasta_fpkm[$i-1][$j];
#		}
#	for(my $i=0;$i<$fpkms;$i++){	# 二次元配列のunshift をどうやるかわからなかったのでむりやり	要らなかった疑惑
#		$fasta_fpkm[0][$i] = $hitcontig_fpkm[$hitnum][$i];)
#	}
#	search_fasta_fpkm();	# 二次元配列のunshiftが面倒だったので全リスト検索しちゃう
#	if($trinities>0){search_fasta_seq();}
#		$textwidget -> delete('1.0','end');
#		for(my $i=0;$i<$fastas;$i++){
#			$textwidget -> insert('end','>'.$fasta_name[$i]."\n");
#			$textwidget -> insert('end',$fasta_seq[$i]."\n");
#		}
#	$textwidget_status = 'available';
#	draw_scat();
	refresh_fasta();
}

#########
# スキャッタープロット上で背景部分がクリックされた時の動作
#########

sub create_by_userclick(){
	# クリックされた背景画像上の座標を調べる
	my @argument = @_;
	my $x = $argument[1];
	my $y = $argument[2];
	$mousex = $plot_canvas -> canvasx($x);
	$mousey = $plot_canvas -> canvasy($y);
	my $z = $scat_zoom;		# なぜかこうしないとバグる。なんでだろう・・
	if($query_length == 0){if($keyword_hits == 0){
		my $sc = int($x/$graphsize);	# 何枚目のグラフか調べる
		my $mergin = 80 * $z;	# 左と下の余白
		my $scale = 100 * $z;	# 10倍目盛のサイズ
		my $xzero = $mergin + $graphsize * $sc;
		my $yzero = $graphsize - $mergin;

		# 一番近い
		my $near=10000;
		my $id;
		for(my $i=0;$i<$contigs;$i++){
			my $x1 = $xzero + $scale * (log($contig_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3);
			my $y1 = $yzero - $scale * (log($contig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3);
			my $dist = sqrt(($x1-$mousex)**2 + ($y1-$mousey)**2);
			if($near > $dist){
					$near = $dist;
					$id   = $i;
			}
		}
		if(defined($id)){
			if($textwidget_refreshflag == 1){input_fasta();}
			$fasta_name[$fastas] = $contig_name[$id];		# 背景クリックが有効なのはBLASTとキーワードがオフのときだけなので後ろに追加でいい
			$fastas++;
#			search_fasta_fpkm();	# fpkm値を探す。$id = 目標contigのIDなので検索可能だが面倒なのでsubを使う
#			$textwidget_refreshflag = 0;
#			search_fasta_seq();		# シーケンスを探す。全コンティグの配列をアサインしてる訳でないのでcontigIDからの検索は無理
#			draw_scat();
			refresh_fasta();
		}
	}}
}

#########
# 散布図を保存するファイルを選ぶ。
#########

sub save_scat(){
	$scat_savefile = $plot_window -> getSaveFile( -filetypes => [['PNG image','.png'],['All files','*']]);
	if(defined($scat_savefile)){
#		draw_scat();
#		$plot_canvas -> postscript(-file => $scat_savefile);
		if($scat_savefile =~ /\.png/){}else{$scat_savefile = $scat_savefile.'.png';}
		open OUT, '>'.$scat_savefile;
		binmode OUT;
		print OUT $saveimage_scat -> png();
		close OUT;
	}else{
		print STDERR 'save canceled'."\n";
	}
	save_settings();
}

#########
# Localblast からTrinity を検索してFastaのリストに追加する。キーワード検索も統合ボタンにする
#########
sub blast2fasta(){
	if($hitcontigs>0){
		@keyword_hit_name = ();
		@keyword_hit_seq = ();
		@keyword_hit_fpkm = ();
		@keyword_hit_length = ();
		$keyword_hits = 0;
		$keyword_error = 'disable';

		if($trinities>0){
			my $count=0;
			for(my $i=0;$i<$hitcontigs;$i++){
					TRI:for(my $j=0;$j<$trinities;$j++){
						if($hitcontig_name[$i] eq $trinity_name[$j]){
							$hitcontig_seq[$i] = $trinity_seq[$j];
							$count++;
							last TRI;
						}
					}
			}
		}else{
			print STDERR 'Require loading Trinity.fasta to get seq'."\n";
		}
			$textwidget2 -> delete('1.0','end');
			for(my $i=0;$i<$hitcontigs;$i++){
				$textwidget2 -> insert('end','>'.$hitcontig_name[$i]."\n");
				$textwidget2 -> insert('end',$hitcontig_seq[$i]."\n");
			}
			print STDERR 'blast2fasta: '.$count.'contigs found'."\n";
			
			open OUT, '>RNAseqViewer_blasthits.fa';		# ヒットした配列をまとめてファイルに出力
			for(my $i=0;$i<$hitcontigs;$i++){
				print OUT '>'.$hitcontig_name[$i]."\n";
				print OUT $hitcontig_seq[$i]."\n";
			}
			close OUT;
			
	}else{

		# キーワード検索結果から配列を検索する。ヒット数が多いと時間がかかる。
		if($query_length == 0){
			if($keyword_hits>0){
				if($trinities>0){
					my $count=0;
					for(my $i=0;$i<$keyword_hits;$i++){
							TRI:for(my $j=0;$j<$trinities;$j++){
								if($keyword_hit_name[$i] eq $trinity_name[$j]){
									$keyword_hit_seq[$i] = $trinity_seq[$j];
									print STDERR '.';
									$count++;
									last TRI;
								}
							}
					}
					print STDERR "\n".'keyword search: sequence found for '.$count.' / '.$keyword_hits." contigs\n";
				}else{
					print STDERR 'Require loading Trinity.fasta to get seq'."\n";
				}
				$textwidget2 -> delete('1.0','end');
				for(my $i=0;$i<$keyword_hits;$i++){
					$textwidget2 -> insert('end','>'.$keyword_hit_name[$i]."\n");
					$textwidget2 -> insert('end',$keyword_hit_seq[$i]."\n");
				}
			}
		}

	}
}


#########
# テキストウィジェットに入力されたデータをメモリーにロードする
#########

sub input_fasta(){
	# まずチェック
	my $text = $textwidget -> get('1.0','end-1c');
	   $text =~ s/ //g;	# まず空白は全て消す
	   $text =~ s/\t//g;	# タブも全て消す
	   $text =~ s/\|/_/g;	# 縦棒は全て_で変換
	my @splitdata = split(/>/,$text);
	my $count = 0;
	my $check = 1;
	CHK:foreach my $entrydata (@splitdata){
		if($count>0){
			my @s = split(/\n/,$entrydata);
			my $name = $s[0];
			   $name =~ s/\r//;		# >で分けたうちの最初の改行までを名前とする
			my $seq = $entrydata;
			   $seq =~ s/\n//g;		# 改行と名前の分を取り除いた残りをシーケンスとして取り扱う
			   $seq =~ s/\r//g;
			   $seq =~ s/$name//;
			if(length($name)>20){				# シーケンスだけ入力されていたりする場合は名前欄が２０文字を超えたりする
				$textwidget_status = 'invalid1';
				$check = 0;
				last CHK;
			}
			if($seq =~ /[0-9]/){				# シーケンスに数字が入っていたら名前が入っちゃってる
				$textwidget_status = 'invalid2';
				$check = 0;
				last CHK;
			}
		}else{
			if(length($entrydata) > 0){			# 最初の行が > で始まってない場合はエラーとする
				$textwidget_status = 'invalid3';
				$check = 0;
				last CHK;
			}
		}
		$count++;
	}

	if($check == 1){
		$fastas=0;
		my $countt=0;
		foreach my $entrydata(@splitdata){
			if($countt>0){
				my @s = split(/\n/,$entrydata);
				my $name = $s[0];
				   $name =~ s/\r//;		# >で分けたうちの最初の改行までを名前とする
				my $seq = $entrydata;
				   $seq =~ s/\n//g;		# 改行と名前の分を取り除いた残りをシーケンスとして取り扱う
				   $seq =~ s/\r//g;
				   $seq =~ s/$name//;
				$fasta_name[$fastas] = $name;
				$fasta_seq[$fastas]  = $seq;
				$fastas++;
			}
			$countt++;
		}
#		$textwidget -> delete('1.0','end');
#		for(my $i=0;$i<$fastas;$i++){
#			$textwidget -> insert('end','>'.$fasta_name[$i]."\n");
#			$textwidget -> insert('end',$fasta_seq[$i]."\n");
#		}
		$textwidget_status = 'available';
		$textwidget_refreshflag = 0;
		refresh_fasta();
	}else{
		#$textwidget_status = 'invalid';
	}
}

#########
# Fastaリストをファイルに保存する
#########

sub save_fasta(){
	#search_fasta_seq();
	if($textwidget_status eq 'available'){
	my $savefile = $control_window -> getSaveFile( -filetypes => [['Fasta','.txt'],['Fasta','.fasta'],['Fasta','.fas'],['All files','*']]);
	if(defined($savefile)){
		open OUT,'>'.$savefile;
			for(my $i=0;$i<$fastas;$i++){
				print OUT '>'.$fasta_name[$i]."\n";
				print OUT $fasta_seq[$i]."\n";
			}
		close OUT;
	}else{
		print STDERR 'save cancelled'."\n";
	}
	}
}
my $entrylist_file = 'contiglist.fasta';		# IDリストのファイル
my $entrylist_file_status = 'not available';# IDリストの読み込み状況

#########
# ファイルからFastaリストをロードする
#########

sub findfile_fasta(){
	my $loadfile = $control_window -> getOpenFile( -filetypes => [['Fasta','.txt'],['Fasta','.fasta'],['Fasta','.fas'],['All files','*']]);
	if(defined($loadfile)){
		$entrylist_file = $loadfile;
		my $seqin = Bio::SeqIO->new(-file => $loadfile, -format => 'fasta');
		while(my $seqobj = $seqin -> next_seq){
			$fasta_name[$fastas] = $seqobj -> id;
			$fasta_name[$fastas] =~ s/\|/_/;
			$fasta_seq[$fastas] = $seqobj -> seq;
			$fastas++;
		}
		$textwidget -> delete('1.0','end');
		for(my $i=0;$i<$fastas;$i++){
			$textwidget -> insert('end','>'.$fasta_name[$i]."\n");
			$textwidget -> insert('end',$fasta_seq[$i]."\n");
		}
		$control_window -> update;
		input_fasta();
		$textwidget_status = 'available';
	}else{
		print STDERR 'load cancelled'."\n";
	}

}

#########
# Fastaリストをクリアする
#########

sub clear_fasta(){
	@fasta_name=();
	@fasta_seq=();
	@fasta_fpkm=();
	$fastas=0;
	$textwidget -> delete('1.0','end');
	$textwidget_status = 'empty';
	draw_scat();
}


#########
# ユーザー指定欄のコンティグ名からfpkm値を検索する
#########

sub search_fasta_fpkm(){
	if($contigs > 0){
		my $count=0;
		for(my $j=0;$j<$fastas;$j++){
			HIT:for(my $i=0;$i<$contigs;$i++){
				if($fasta_name[$j] eq $contig_name[$i]){
					for(my $k=0;$k<$fpkms;$k++){
						$fasta_fpkm[$j][$k] = $contig_fpkm[$i][$k];
					}
					$count++;
					last HIT;
				}
			}
		}
		print STDERR 'search_fasta_fpkm: '.$count.' / '.$hitcontigs."\n";
	}else{
		print STDERR 'RSEM data not found'."\n";
	}
}

#########
# ユーザー指定欄のコンティグ名からTrinity を検索してFastaのリストに追加する
#########
sub search_fasta_seq(){
#	if($textwidget_refreshflag == 1){input_fasta();}
#	if($textwidget_status eq 'available'){
	if($fastas>0){
		if($trinities>0){
			my $count=0;
			for(my $i=0;$i<$fastas;$i++){
				if(length($fasta_seq[$i])==0){
					TRI:for(my $j=0;$j<$trinities;$j++){
						if($fasta_name[$i] eq $trinity_name[$j]){
							$fasta_seq[$i] = $trinity_seq[$j];
							$count++;
							last TRI;
						}
					}
				}
			}
			print STDERR 'search_fasta_seq: '.$count.'contigs found'."\n";
		}else{
			print STDERR 'Trinity.fasta have not loaded; necessary to find contig sequence'."\n";
		}
#			$textwidget -> delete('1.0','end');
#			for(my $i=0;$i<$fastas;$i++){
#				$textwidget -> insert('end','>'.$fasta_name[$i]."\n");
#				$textwidget -> insert('end',$fasta_seq[$i]."\n");
#			}
	}
#	}
}

#########
# ユーザー指定欄をリフレッシュ
#########
sub refresh_fasta(){
		$textwidget -> delete('1.0','end');
		search_fasta_fpkm();
		search_fasta_seq();
		for(my $i=0;$i<$fastas;$i++){
			$textwidget -> insert('end','>'.$fasta_name[$i]."\n");
			$textwidget -> insert('end',$fasta_seq[$i]."\n");
		}
		$textwidget_status = 'available';
		draw_scat();

}
