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
    public interface IQtlAnalysisScenarioSettingValue : IQtlSeqTargetPolicySettingValue, INoQtlDistributionOption, ISlidingWindowAnalysisOption, IGraphOptions
    {
        /// <summary>
        /// 出力ディレクトリを取得する。
        /// </summary>
        string OutputDir { get; }

        /// <summary>
        /// 表示するAnnotation Imapctを指定する。
        /// </summary>
        string DisplayAnnotationImpacts { get; }

        /// <summary>
        /// スレッド数を指定する。
        /// </summary>
        int ThreadNumber { get; }
    }
}
