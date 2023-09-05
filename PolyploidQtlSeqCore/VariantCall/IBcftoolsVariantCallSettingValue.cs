namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// Bcftools変異検出オプションインターフェイス
    /// </summary>
    public interface IBcftoolsVariantCallSettingValue
    {
        /// <summary>
        /// リファレンスシークエンスファイルを取得する。
        /// </summary>
        string ReferenceSequence { get; }

        /// <summary>
        /// MQ最低値を取得する。
        /// </summary>
        int MinMq { get; }

        /// <summary>
        /// BQ最低値を取得する。
        /// </summary>
        int MinBq { get; }

        /// <summary>
        /// Adjust MQを取得する。
        /// </summary>
        int AdjustMq { get; }

        /// <summary>
        /// 出力ディレクトリを取得する。
        /// </summary>
        string OutputDir { get; }

        /// <summary>
        /// スレッド数を取得する。
        /// </summary>
        int ThreadNumber { get; }
    }
}
