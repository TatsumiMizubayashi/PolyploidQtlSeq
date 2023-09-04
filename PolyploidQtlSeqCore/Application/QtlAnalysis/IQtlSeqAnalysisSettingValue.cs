namespace PolyploidQtlSeqCore.Application.QtlAnalysis
{
    /// <summary>
    /// QTL-seq解析設定値インターフェース
    /// </summary>
    public interface IQtlSeqAnalysisSettingValue
    {
        /// <summary>
        /// InputVCFファイルPathを取得、又は設定する。
        /// </summary>
        string InputVcf { get; set; }

        /// <summary>
        /// 出力ディレクトリを取得、又は設定する。
        /// </summary>
        string OutputDir { get; set; }

        /// <summary>
        /// 表示するAnnotation Imapctを取得、又は設定する。
        /// </summary>
        string DisplayAnnotationImpacts { get; set; }

        /// <summary>
        /// スレッド数を取得、又は設定する。
        /// </summary>
        int ThreadNumber { get; set; }

        /// <summary>
        /// P1 MostAlleleRateのしきい値を取得、又は設定する。
        /// </summary>
        double Parent1MostAlleleRateThreshold { get; set; }

        /// <summary>
        /// P2 SNP-index範囲を取得、又は設定する。
        /// </summary>
        string Parent2SnpIndexRange { get; set; }

        /// <summary>
        /// 最低Depthのしきい値を取得、又は設定する。
        /// </summary>
        int MinimumDepthThreshold { get; set; }

        /// <summary>
        /// 最大Bulk SNP-indexのしきい値を取得、又は設定する。
        /// </summary>
        double MaxBulkSnpIndexThreshold { get; set; }

        /// <summary>
        /// 倍数性の値を取得、又は設定する。
        /// </summary>
        int Ploidy { get; set; }

        /// <summary>
        /// P2のplex数を取得、又は設定する。
        /// </summary>
        int Parent2PlexNumber { get; set; }

        /// <summary>
        /// Bulk1の個体数を取得、又は設定する。
        /// </summary>
        int Bulk1Number { get; set; }

        /// <summary>
        /// Bulk2の個体数を取得、又は設定する。
        /// </summary>
        int Bulk2Number { get; set; }

        /// <summary>
        /// 分布作成時の試行回数を取得、又は設定する。
        /// </summary>
        int ReplicatesNumber { get; set; }

        /// <summary>
        /// Windowサイズ(kbp)を取得、又は設定する。
        /// </summary>
        int WindowSize { get; set; }

        /// <summary>
        /// Window移動量(kbp)を取得、又は設定する。
        /// </summary>
        int StepSize { get; set; }

        /// <summary>
        /// グラフ画像の幅(Pixel)を取得、又は設定する。
        /// </summary>
        int FigureWidth { get; set; }

        /// <summary>
        /// グラフ画像の高さ(Pixel)を取得、又は設定する。
        /// </summary>
        int FigureHeight { get; set; }

        /// <summary>
        /// X軸の目盛り間隔(MB)を取得、又は設定する。
        /// </summary>
        int XAxisMajorStep { get; set; }
    }
}
