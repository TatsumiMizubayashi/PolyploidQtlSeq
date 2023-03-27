using PolyploidQtlSeqCore.QtlAnalysis.IO;

namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// 最大スコアのSlidingWindowQTl情報
    /// </summary>
    internal class MaxScoreSlidingWindowQtl
    {
        /// <summary>
        /// 最大スコアのSlidingWindowQTl情報を作成する。
        /// </summary>
        /// <param name="window">最大スコアWindow</param>
        public MaxScoreSlidingWindowQtl(Window window)
        {
            Score = window.AverageScore.Value;
            PValue = window.AveragePValue.Value;
            IsP95Qtl = window.P95Qtl.IsQtl;
            IsP99Qtl = window.P99Qtl.IsQtl;
        }

        /// <summary>
        /// 最大スコアを取得する。
        /// </summary>
        public double Score { get; }

        /// <summary>
        /// 最低PValueを取得する。
        /// </summary>
        public double PValue { get; }

        /// <summary>
        /// P95でQTLがあるかどうかを取得する。
        /// </summary>
        public bool IsP95Qtl { get; }

        /// <summary>
        /// P99でQTLがあるかどうかを取得する。
        /// </summary>
        public bool IsP99Qtl { get; }

        /// <summary>
        /// P95QTLの印を取得する。
        /// </summary>
        /// <returns>P95QTLの印</returns>
        public string P95QtlMark()
        {
            return IsP95Qtl
                ? SnpIndexFileCreator.QTL
                : "";
        }

        /// <summary>
        /// P99QTLの印を取得する。
        /// </summary>
        /// <returns>P99QTLの印</returns>
        public string P99QtlMark()
        {
            return IsP99Qtl
                ? SnpIndexFileCreator.QTL
                : "";
        }
    }
}
