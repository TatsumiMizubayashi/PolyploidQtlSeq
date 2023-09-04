using PolyploidQtlSeqCore.QtlAnalysis.Distribution;
using PolyploidQtlSeqCore.QtlAnalysis.OxyGraph;
using PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter;
using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// QTL解析シナリオオプション インターフェイス
    /// </summary>
    [Obsolete("internalにする予定")]
    public interface IQtlAnalysisScenarioSettingValue : IQtlSeqTargetPolicySettingValue, INoQtlDistributionSettingValue, 
        ISlidingWindowAnalysisSettingValue, IGraphSettingValue
    {
        /// <summary>
        /// 出力ディレクトリを取得する。
        /// </summary>
        string OutputDir { get; }

        /// <summary>
        /// 表示するAnnotation Imapctを取得する。
        /// </summary>
        string DisplayAnnotationImpacts { get; }

        /// <summary>
        /// スレッド数を取得する。
        /// </summary>
        int ThreadNumber { get; }
    }
}
