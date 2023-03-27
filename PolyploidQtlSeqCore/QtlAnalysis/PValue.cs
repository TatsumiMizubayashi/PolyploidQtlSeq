namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// PValue
    /// </summary>
    internal class PValue
    {
        /// <summary>
        /// PValueを作成する。
        /// </summary>
        /// <param name="value">P値</param>
        public PValue(double value)
        {
            if (value <= 0 || value > 1) throw new ArgumentException(null, nameof(value));

            Value = value;
        }

        /// <summary>
        /// PValueを取得する。
        /// </summary>
        public double Value { get; }

        /// <summary>
        /// スコアに変換する。
        /// </summary>
        /// <returns>スコア</returns>
        public Score ToScore()
        {
            var score = -1 * Math.Log10(Value);

            return new Score(score);
        }
    }
}
