using PolyploidQtlSeqCore.QtlAnalysis;

namespace PolyploidQtlSeqCore.Application.QtlAnalysis
{
    /// <summary>
    /// QTL解析コマンドオプション インターフェイス
    /// </summary>
    public interface IQtlAnalysisCommandOptions : IQtlAnalysisScenarioOptions
    {
        /// <summary>
        /// InputVCFファイルPathを取得する。
        /// </summary>
        string InputVcf { get; }

        /// <summary>
        /// パラメータファイルPathを取得する。
        /// </summary>
        string ParameterFile { get; }
    }
}
