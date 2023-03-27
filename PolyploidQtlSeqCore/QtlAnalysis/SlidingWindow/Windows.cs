namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// ウインドウコレクション
    /// </summary>
    internal class Windows
    {
        private readonly Window[] _windows;

        /// <summary>
        /// ウインドウコレクションを作成する。
        /// </summary>
        /// <param name="windows">ウインドウ配列</param>
        public Windows(Window[] windows)
        {
            _windows = windows;
        }

        /// <summary>
        /// 最大スコアWindowQTL情報に変換する。
        /// </summary>
        /// <returns>最大スコアWindowQTL情報</returns>
        public MaxScoreSlidingWindowQtl ToMaxScoreSlidingWindowQtl()
        {
            var maxScoreWindow = _windows.MaxBy(x => x.AverageScore.Value);
            if (maxScoreWindow == null) throw new ArgumentNullException(nameof(maxScoreWindow));

            return new MaxScoreSlidingWindowQtl(maxScoreWindow);
        }
    }
}
