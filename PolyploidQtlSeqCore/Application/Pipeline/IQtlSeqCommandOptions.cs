using PolyploidQtlSeqCore.Application.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;
using PolyploidQtlSeqCore.QtlAnalysis.Mapping;
using PolyploidQtlSeqCore.QtlAnalysis.VariantCall;

namespace PolyploidQtlSeqCore.Application.Pipeline
{
    /// <summary>
    /// QTL-seq解析コマンドオプションインターフェイス
    /// </summary>
    public interface IQtlSeqCommandOptions : IMappingOptions, IAnalysisChrOptions, IBcftoolsVariantCallOption,
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
