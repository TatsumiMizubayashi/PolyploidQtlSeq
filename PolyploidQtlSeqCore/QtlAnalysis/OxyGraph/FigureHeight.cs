namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ画像の高さ
    /// </summary>
    internal class FigureHeight
    {
        /// <summary>
        /// グラフ高さの最小値
        /// </summary>
        private const int MINIMUM = 100;

        /// <summary>
        /// グラフ高さの最大値
        /// </summary>
        private const int MAXIMUM = 2000;


        /// <summary>
        /// グラフ画像の高さ(pixel)を作成する。
        /// </summary>
        /// <param name="height">高さ(pixel)</param>
        public FigureHeight(int height)
        {
            if (height   < MINIMUM || height > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(height));

            Value = height;
        }

        /// <summary>
        /// グラフ画像の高さ(Pixel)を取得する。
        /// </summary>
        internal int Value { get; }
    }
}
