namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEff オプションインターフェイス
    /// </summary>
    public interface ISnpEffOption
    {
        /// <summary>
        /// SnpEffに指定するjava Max Heapサイズ(GB)を取得する。
        /// </summary>
        int SnpEffMaxHeap { get; }

        /// <summary>
        /// SnpEff.configファイルのPathを取得する。
        /// </summary>
        string SnpEffConfigFile { get; }

        /// <summary>
        /// SnpEffで使用するデータベース名を取得する。
        /// </summary>
        string SnpEffDatabaseName { get; }
    }
}
