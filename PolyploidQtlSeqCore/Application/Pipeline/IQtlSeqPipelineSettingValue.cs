using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;
using PolyploidQtlSeqCore.VariantCall;

namespace PolyploidQtlSeqCore.Application.Pipeline
{
    /// <summary>
    /// QTL-seq解析パイプライン設定値 インターフェース
    /// </summary>
    public interface IQtlSeqPipelineSettingValue : IMappingOptions, IAnalysisChrOptions, IBcftoolsVariantCallOption,
        ISnpEffOption, IQtlAnalysisScenarioOptions
    {
        /// <summary>
        /// リファレンスシークエンスファイルを取得する。
        /// </summary>
        string ReferenceSequence { get; }

        /// <summary>
        /// パラメータファイルを取得する。
        /// </summary>
        string ParameterFile { get; }
    }
}
