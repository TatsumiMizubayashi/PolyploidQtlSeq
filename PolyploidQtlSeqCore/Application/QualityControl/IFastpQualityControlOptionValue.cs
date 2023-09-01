namespace PolyploidQtlSeqCore.Application.QualityControl
{
    /// <summary>
    /// Fastp QCオプション値 インターフェース
    /// </summary>
    public interface IFastpQualityControlOptionValue
    {
        /// <summary>
        /// 入力Fastqファイル置き場のPathを設定、または取得する。
        /// </summary>
        string InputDir { get; set; }

        /// <summary>
        /// 出力ディレクトリのPathを設定、または取得する。
        /// </summary>
        string OutputDir { get; set; }


        /// <summary>
        /// トリミング後リード長最低値を設定、または取得する。
        /// </summary>
        int ReadLengthRequired { get; set; }

        /// <summary>
        /// N塩基数の上限を設定、または取得する。
        /// </summary>
        int NBaseLimit { get; set; }

        /// <summary>
        /// 塩基クオリティのしきい値を設定、または取得する。
        /// </summary>
        int BaseQuality { get; set; }

        /// <summary>
        /// 3'末端トリム時の平均クオリティしきい値を設定、または取得する。
        /// </summary>
        int CutTailMeanQuality { get; set; }

        /// <summary>
        /// 3'末端トリム時のウインドウサイズを設定、または取得する。
        /// </summary>
        int CutTailWindowSize { get; set; }

        /// <summary>
        /// 使用するスレッド数を設定、または取得する。
        /// </summary>
        int ThreadNumber { get; set; }
    }
}
