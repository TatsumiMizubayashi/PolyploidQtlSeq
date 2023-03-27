using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// Window Score Y軸設定クリエーター
    /// </summary>
    internal static class WindowScoreYAxisConfigCreator
    {
        private static readonly double _padding = 0.1;
        private static readonly double _step = 1;

        /// <summary>
        /// Window ScoreのY軸設定を作成する。
        /// </summary>
        /// <param name="windows"></param>
        /// <returns>Window ScoreのY軸設定</returns>
        public static YAxisConfig Create(Window[] windows)
        {
            var maxScore = windows.Max(x => x.AverageScore.Value);

            var max = maxScore + _padding;

            return new YAxisConfig(0, max, _step);
        }
    }
}
