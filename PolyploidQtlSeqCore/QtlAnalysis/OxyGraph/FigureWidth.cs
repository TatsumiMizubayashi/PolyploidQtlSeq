namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ画像の幅
    /// </summary>
    internal class FigureWidth
    {
        /// <summary>
        /// グラフ幅の最小値
        /// </summary>
        private const int MINIMUM = 300;

        /// <summary>
        /// グラフ幅の最大値
        /// </summary>
        private const int MAXIMUM = 5000;

        /// <summary>
        /// グラフ画像の幅を作成する。
        /// </summary>
        /// <param name="width">幅(pixel)</param>
        public FigureWidth(int width)
        {
            if (width < MINIMUM || width > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(width));

            Value = width;
        }

        /// <summary>
        /// グラフ画像の幅(Pixel)を取得する。
        /// </summary>
        internal int Value { get; }

    }
}
