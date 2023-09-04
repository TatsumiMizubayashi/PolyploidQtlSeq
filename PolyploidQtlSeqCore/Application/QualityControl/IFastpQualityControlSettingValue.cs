namespace PolyploidQtlSeqCore.Application.QualityControl
{
    /// <summary>
    /// Fastp QCオプション値 インターフェース
    /// </summary>
    public interface IFastpQualityControlSettingValue
    {
        /// <summary>
        /// 入力Fastqファイル置き場のPathを取得する。
        /// </summary>
        string InputDir { get; }

        /// <summary>
        /// 出力ディレクトリのPathを取得する。
        /// </summary>
        string OutputDir { get; }


        /// <summary>
        /// トリミング後リード長最低値を取得する。
        /// </summary>
        int ReadLengthRequired { get; }

        /// <summary>
        /// N塩基数の上限を取得する。
        /// </summary>
        int NBaseLimit { get; }

        /// <summary>
        /// 塩基クオリティのしきい値を取得する。
        /// </summary>
        int BaseQuality { get; }

        /// <summary>
        /// 3'末端トリム時の平均クオリティしきい値を取得する。
        /// </summary>
        int CutTailMeanQuality { get; }

        /// <summary>
        /// 3'末端トリム時のウインドウサイズを取得する。
        /// </summary>
        int CutTailWindowSize { get; }

        /// <summary>
        /// 使用するスレッド数を取得する。
        /// </summary>
        int ThreadNumber { get; }
    }
}
