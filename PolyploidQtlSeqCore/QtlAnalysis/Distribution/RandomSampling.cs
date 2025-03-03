namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// ランダムサンプリング
    /// </summary>
    internal class RandomSampling
    {
        private readonly Random _random;

        /// <summary>
        /// ランダムサンプリングインスタンスを作成する。
        /// </summary>
        public RandomSampling()
        {
            _random = new Random();
        }

        /// <summary>
        /// 重複を許す無作為抽出を行う。
        /// </summary>
        /// <param name="values">抽出対象値の配列</param>
        /// <param name="number">取り出す数</param>
        /// <returns>抽出値の配列</returns>
        public int[] Choice(int[] values, int number)
        {
            var indexes = new int[number];
            for (var i = 0; i < number; i++)
            {
                indexes[i] = _random.Next(0, values.Length);
            }

            return [.. indexes.Select(i => values[i])];
        }

        /// <summary>
        /// 重みを付けて重複を許す無作為抽出を行う。
        /// </summary>
        /// <param name="values">抽出対象値の配列</param>
        /// <param name="weights">重み</param>
        /// <param name="number">取り出す数</param>
        /// <returns>抽出値の配列</returns>
        public int[] Choice(int[] values, double[] weights, int number)
        {
            var weightTotal = weights.Sum();
            if (weightTotal == 0.0) throw new ArgumentException($"The total value of {nameof(weights)} is 0.");

            var items = values.Zip(weights, (v, w) => new { Value = v, Weight = w / weightTotal }).ToArray();
            var choiceValues = new int[number];
            for (var i = 0; i < number; i++)
            {
                var randValue = _random.NextDouble();
                for (var j = 0; j < items.Length; j++)
                {
                    var item = items[j];
                    if (item.Weight == 0) continue;
                    if (item.Weight > randValue)
                    {
                        choiceValues[i] = item.Value;
                        break;
                    }

                    randValue -= item.Weight;
                }
            }

            return choiceValues;
        }
    }
}
