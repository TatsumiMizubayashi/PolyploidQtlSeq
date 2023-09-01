using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;
using PolyploidQtlSeqCore.VariantCall;

namespace PolyploidQtlSeqCore.Application.Pipeline
{
    /// <summary>
    /// QTL-seq解析パイプライン設定値 インターフェース
    /// </summary>
    public interface IQtlSeqPipelineSettingValue : IMappingSettingValue, IMappingSampleSettingValue, IAnalysisChrSettingValue, 
        IBcftoolsVariantCallSettingValue, ISnpEffSettingValue, IQtlAnalysisScenarioSettingValue
    {
        /// <summary>
        /// パラメータファイルを取得する。
        /// </summary>
        string ParameterFile { get; }
    }
}
