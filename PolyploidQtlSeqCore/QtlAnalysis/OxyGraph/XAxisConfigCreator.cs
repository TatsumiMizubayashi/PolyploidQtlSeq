using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// ゲノム位置 X軸設定クリエーター
    /// </summary>
    internal static class XAxisConfigCreator
    {
        /// <summary>
        /// X軸設定を作成する。
        /// </summary>
        /// <param name="windows">全Window</param>
        /// <param name="majorStep">X軸MajorStep(MB)</param>
        /// <returns>X軸設定</returns>
        public static XAxisConfig Create(Window[] windows, XAxisMajorStep majorStep)
        {
            var maxBp = windows
                .Where(x => x.VariantCount.Count > 0)
                .Max(x => x.GenomePosition.End);

            var mbp = maxBp.ToMbp();
            var celingMaxMbp = Ceiling(mbp);

            return new XAxisConfig(0, celingMaxMbp, majorStep);
        }

        /// <summary>
        /// Mbpの値を小数点第2位で切り上げる。
        /// </summary>
        /// <param name="mbp">Mbp</param>
        /// <returns>小数点第２位で切り上げたMbp</returns>
        private static double Ceiling(double mbp)
        {
            var tmp = mbp * 10;
            var celingTmp = Math.Ceiling(tmp);

            return celingTmp / 10;
        }
    }
}
