using PolyploidQtlSeqCore.MappingAndVariantCall;
using PolyploidQtlSeqCore.QtlAnalysis;

namespace PolyploidQtlSeqCore.Application.Pipeline
{
    /// <summary>
    /// QTL-seq解析パイプライン設定値 インターフェース
    /// </summary>
    public interface IQtlSeqPipelineSettingValue : IVariantCallPipelineSettingValue, IQtlAnalysisScenarioSettingValue
    {
        /// <summary>
        /// パラメータファイルを取得する。
        /// </summary>
        [Obsolete("削除予定")]
        string ParameterFile { get; }
    }
}
