namespace PolyploidQtlSeqCore.QtlAnalysis.Mapping
{
    /// <summary>
    /// Mappingオプションインターフェイス
    /// </summary>
    public interface IMappingOptions
    {
        /// <summary>
        /// 親1ディレクトリを取得する。
        /// </summary>
        string Parent1Dir { get; }

        /// <summary>
        /// 親2ディレクトリを取得する。
        /// </summary>
        string Parent2Dir { get; }

        /// <summary>
        /// Bulk1ディレクトリを取得する。
        /// </summary>
        string Bulk1Dir { get; }

        /// <summary>
        /// Bulk2ディレクトリを取得する。
        /// </summary>
        string Bulk2Dir { get; }
    }
}
