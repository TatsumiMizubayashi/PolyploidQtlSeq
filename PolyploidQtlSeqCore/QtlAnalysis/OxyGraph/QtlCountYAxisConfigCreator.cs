using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// QTL数 Y軸設定クリエーター
    /// </summary>
    internal static class QtlCountYAxisConfigCreator
    {
        private static readonly int _padding = 5;

        /// <summary>
        /// QTL数 Y軸設定を作成する。
        /// </summary>
        /// <param name="windows">windows</param>
        /// <returns>Y軸設定</returns>
        public static YAxisConfig Create(Window[] windows)
        {
            var maxP95PlusCount = windows.Max(x => x.P95Qtl.PlusDeltaSnpIndexQtlVariantCount);
            var maxP95MinusCount = windows.Max(x => x.P95Qtl.MinusDeltaSnpIndexQtlVariantCount);
            var step = GetStep(maxP95PlusCount, maxP95MinusCount);

            // QTLグラフではΔSNP-indexが負の値はカウント*-1として表現する
            var min = -1 * (maxP95MinusCount + _padding);
            var max = maxP95PlusCount + _padding;

            return new YAxisConfig(min, max, step);
        }

        /// <summary>
        /// Y軸Setpを取得する。
        /// </summary>
        /// <param name="plusCount">Plus QTL数</param>
        /// <param name="minusCount">Minus QTL数</param>
        /// <returns>Step</returns>
        private static int GetStep(int plusCount, int minusCount)
        {
            var length = plusCount + minusCount;

            if (length <= 150) return 25;
            if (length <= 300) return 50;

            return 100;
        }
    }
}
