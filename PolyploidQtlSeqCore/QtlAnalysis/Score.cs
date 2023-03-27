namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// スコア
    /// </summary>
    internal class Score
    {
        /// <summary>
        /// スコアを作成する。
        /// </summary>
        /// <param name="value">スコア</param>
        public Score(double value)
        {
            if (value < 0) throw new ArgumentException(null, nameof(value));

            Value = value;
        }

        /// <summary>
        /// スコアを取得する。
        /// </summary>
        public double Value { get; }

        /// <summary>
        /// PValueに変換する。
        /// </summary>
        /// <returns>PValue</returns>
        public PValue ToPValue()
        {
            var p = Math.Pow(10, -1 * Value);
            return new PValue(p);
        }
    }
}
