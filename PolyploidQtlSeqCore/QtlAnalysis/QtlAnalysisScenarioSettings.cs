using PolyploidQtlSeqCore.QtlAnalysis.Distribution;
using PolyploidQtlSeqCore.QtlAnalysis.IO;
using PolyploidQtlSeqCore.QtlAnalysis.OxyGraph;
using PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter;
using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// QTL解析シナリオ設定
    /// </summary>
    internal class QtlAnalysisScenarioSettings
    {
        /// <summary>
        /// QTL解析しなりを設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値</param>
        public QtlAnalysisScenarioSettings(IQtlAnalysisScenarioSettingValue settingValue)
        {
            OutputDir = new OutputDirectory(settingValue.OutputDir);
            DisplayAnnotationImpacts = new DisplayAnnotationImpacts(settingValue.DisplayAnnotationImpacts);
            ThreadNumber = new ThreadNumber(settingValue.ThreadNumber);

            QtlSeqTargetPolicySettings = new QtlSeqTargetPolicySettings(settingValue);
            NoQtlDistributionSettings = new NoQtlDistributionSettings(settingValue);
            SlidingWindowAnalysisSettings = new SlidingWindowAnalysisSettings(settingValue);
            GraphSettings = new GraphSettings(settingValue);
        }

        /// <summary>
        /// 出力ディレクトリ
        /// </summary>
        public OutputDirectory OutputDir { get; }

        /// <summary>
        /// 表示するAnnotation Impact
        /// </summary>
        public DisplayAnnotationImpacts DisplayAnnotationImpacts { get; }

        /// <summary>
        /// スレッド数
        /// </summary>
        public ThreadNumber ThreadNumber { get; }

        /// <summary>
        /// QTL-seq対象ポリシー設定
        /// </summary>
        public QtlSeqTargetPolicySettings QtlSeqTargetPolicySettings { get; }

        /// <summary>
        /// QTLなし分布作成設定
        /// </summary>
        public NoQtlDistributionSettings NoQtlDistributionSettings { get; }

        /// <summary>
        /// SlidingWindow解析設定
        /// </summary>
        public SlidingWindowAnalysisSettings SlidingWindowAnalysisSettings { get; }

        /// <summary>
        /// グラフ設定
        /// </summary>
        public GraphSettings GraphSettings { get; }
    }
}
